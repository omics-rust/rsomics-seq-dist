//! Pairwise nucleotide distances matching scikit-bio 0.7.2
//! `skbio.sequence.distance.{hamming, jc69, k2p}`.

use std::str::FromStr;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Metric {
    Hamming,
    Jc69,
    K2p,
}

impl FromStr for Metric {
    type Err = RsomicsError;
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "hamming" => Ok(Self::Hamming),
            "jc69" => Ok(Self::Jc69),
            "k2p" => Ok(Self::K2p),
            other => Err(RsomicsError::InvalidInput(format!(
                "unknown metric {other:?}; expected hamming, jc69, or k2p"
            ))),
        }
    }
}

/// Canonical nucleobase code for K2P bitwise transition/transversion logic:
/// A=0, C=1, G=2, T/U=3. A purine↔purine or pyrimidine↔pyrimidine change XORs
/// to 2 (transition); any other change has an odd XOR (transversion). Code 4 is
/// any non-canonical character (gap, N, ambiguous, lowercase) and is excluded
/// from the valid-site denominator — mirroring scikit-bio's `alphabet="canonical"`
/// pairwise deletion.
const NON_CANONICAL: u8 = 4;

#[inline]
fn nucl_code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => NON_CANONICAL,
    }
}

/// Pre-encode a sequence to canonical codes (0-3 for A/C/G/T-or-U, 4 otherwise)
/// once, so the pairwise hot loop is a branchless scan over `u8` codes the
/// compiler can vectorize, instead of re-classifying every byte per pair.
#[must_use]
pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| nucl_code(b)).collect()
}

/// Per-pair canonical-site tallies over two pre-encoded code slices: number of
/// sites where both are canonical (`valid`), of those the transitions (`ts`,
/// code XOR == 2) and transversions (`tv`, odd XOR). A single branchless pass
/// the optimizer auto-vectorizes; the JC69 differing-site count is `ts + tv`.
#[inline]
#[must_use]
pub fn tally_codes(a: &[u8], b: &[u8]) -> (u64, u64, u64) {
    debug_assert_eq!(a.len(), b.len());
    let mut valid = 0u64;
    let mut ts = 0u64;
    let mut tv = 0u64;
    for (&cx, &cy) in a.iter().zip(b) {
        let both = u64::from(cx != NON_CANONICAL && cy != NON_CANONICAL);
        valid += both;
        let sub = cx ^ cy;
        ts += both & u64::from(sub == 2);
        tv += both & u64::from(sub & 1 == 1);
    }
    (valid, ts, tv)
}

/// JC69 from a pre-tallied `(valid, differing)` count.
#[must_use]
pub fn jc69_from_counts(valid: u64, differing: u64) -> f64 {
    if valid == 0 {
        return f64::NAN;
    }
    let p = differing as f64 / valid as f64;
    if p >= 0.75 {
        return f64::NAN;
    }
    let mut d = p / -0.75;
    d += 1.0;
    d = d.ln();
    d *= -0.75;
    d + 0.0
}

/// K2P from pre-tallied `(valid, transitions, transversions)` counts.
#[must_use]
pub fn k2p_from_counts(valid: u64, ts: u64, tv: u64) -> f64 {
    if valid == 0 {
        return f64::NAN;
    }
    let l = valid as f64;
    let p = ts as f64 / l;
    let q = tv as f64 / l;
    let a1 = 1.0 - 2.0 * p - q;
    let a2 = 1.0 - 2.0 * q;
    if a1 <= 0.0 || a2 <= 0.0 {
        return f64::NAN;
    }
    let d = -0.5 * a1.ln() - 0.25 * a2.ln();
    d + 0.0
}

/// Hamming proportion over all aligned positions, scikit-bio
/// `hamming(proportion=True)`: every character — including gaps and ambiguity
/// codes — is compared literally and the denominator is the full alignment
/// length. Empty input → NaN.
#[must_use]
pub fn hamming(a: &[u8], b: &[u8]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let npos = a.len();
    if npos == 0 {
        return f64::NAN;
    }
    let mismatches = a.iter().zip(b).filter(|(x, y)| x != y).count();
    mismatches as f64 / npos as f64
}

/// p-distance over canonical-only sites (used by JC69): positions where either
/// sequence is non-canonical are dropped, then `mismatches / valid_sites`.
/// No valid site → NaN.
fn p_distance_canonical(a: &[u8], b: &[u8]) -> f64 {
    let mut valid = 0u64;
    let mut diff = 0u64;
    for (&x, &y) in a.iter().zip(b) {
        let cx = nucl_code(x);
        let cy = nucl_code(y);
        if cx == NON_CANONICAL || cy == NON_CANONICAL {
            continue;
        }
        valid += 1;
        if cx != cy {
            diff += 1;
        }
    }
    if valid == 0 {
        return f64::NAN;
    }
    diff as f64 / valid as f64
}

/// JC69 distance, scikit-bio `jc69`: `D = -3/4 · ln(1 - 4/3·p)` over the
/// canonical p-distance. `p == 0 → 0`; `p ≥ 3/4` is over-saturated → NaN
/// (scikit-bio sets `dists[dists >= 0.75] = nan` before the log).
#[must_use]
pub fn jc69(a: &[u8], b: &[u8]) -> f64 {
    let p = p_distance_canonical(a, b);
    if p.is_nan() {
        return f64::NAN;
    }
    if p >= 0.75 {
        return f64::NAN;
    }
    // scikit-bio _p_correct: dists/=-0.75; +=1; ln; *=-0.75; += 0.0 (clear -0.0).
    let mut d = p / -0.75;
    d += 1.0;
    d = d.ln();
    d *= -0.75;
    d + 0.0
}

/// K2P distance, scikit-bio `k2p`:
/// `D = -1/2·ln(1 - 2P - Q) - 1/4·ln(1 - 2Q)` where `P` and `Q` are the
/// transition and transversion proportions over canonical-only sites.
/// Either log argument `≤ 0` (over-saturation) → NaN; identical sequences → 0.
#[must_use]
pub fn k2p(a: &[u8], b: &[u8]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut valid = 0u64;
    let mut ts = 0u64;
    let mut tv = 0u64;
    for (&x, &y) in a.iter().zip(b) {
        let cx = nucl_code(x);
        let cy = nucl_code(y);
        if cx == NON_CANONICAL || cy == NON_CANONICAL {
            continue;
        }
        valid += 1;
        let sub = cx ^ cy;
        if sub == 2 {
            ts += 1;
        } else if sub & 1 == 1 {
            tv += 1;
        }
    }
    if valid == 0 {
        return f64::NAN;
    }
    let l = valid as f64;
    let p = ts as f64 / l;
    let q = tv as f64 / l;

    let a1 = 1.0 - 2.0 * p - q;
    let a2 = 1.0 - 2.0 * q;
    if a1 <= 0.0 || a2 <= 0.0 {
        return f64::NAN;
    }
    let d = -0.5 * a1.ln() - 0.25 * a2.ln();
    d + 0.0
}

#[must_use]
pub fn distance(metric: Metric, a: &[u8], b: &[u8]) -> f64 {
    match metric {
        Metric::Hamming => hamming(a, b),
        Metric::Jc69 => jc69(a, b),
        Metric::K2p => k2p(a, b),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-15;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= EPS, "{a} != {b} (diff {})", (a - b).abs());
    }

    #[test]
    fn hamming_proportion_basic() {
        // scikit-bio doctest: AGGGTA vs CGTTTA -> 0.5
        close(hamming(b"AGGGTA", b"CGTTTA"), 0.5);
        close(hamming(b"ACGT", b"ACGT"), 0.0);
    }

    #[test]
    fn hamming_counts_gaps_literally() {
        // gap '-' vs base is a literal mismatch for hamming; denominator is full length.
        close(hamming(b"ACGT", b"AC-T"), 0.25);
    }

    #[test]
    fn hamming_empty_is_nan() {
        assert!(hamming(b"", b"").is_nan());
    }

    #[test]
    fn jc69_known_value() {
        // s1=ACGTACGTACGT vs s3=GCATACGTACGT, p=2/12 -> 0.18848582121067953
        close(
            jc69(b"ACGTACGTACGT", b"GCATACGTACGT"),
            0.188_485_821_210_679_53,
        );
    }

    #[test]
    fn jc69_identical_is_zero() {
        let d = jc69(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(d, 0.0);
        assert!(d.is_sign_positive(), "must be +0.0 not -0.0");
    }

    #[test]
    fn jc69_saturation_is_nan() {
        // p == 1.0 >= 0.75 -> NaN, no panic on the ln of a negative argument.
        assert!(jc69(b"AAAAAAAA", b"GGGGGGGG").is_nan());
    }

    #[test]
    fn jc69_gap_exclusion_changes_result() {
        // s1 vs s4 (gaps at two positions): canonical filter drops both gap
        // columns; the remaining 10 sites are identical -> distance 0, even
        // though hamming over all 12 chars is non-zero.
        close(jc69(b"ACGTACGTACGT", b"ACGTACG-AC.T"), 0.0);
        assert!(hamming(b"ACGTACGTACGT", b"ACGTACG-AC.T") > 0.0);
    }

    #[test]
    fn k2p_known_value_transitions_only() {
        // s1 vs s3: 2 transitions, 0 transversions -> 0.20273255405408214
        close(
            k2p(b"ACGTACGTACGT", b"GCATACGTACGT"),
            0.202_732_554_054_082_14,
        );
    }

    #[test]
    fn k2p_mixed_transitions_and_transversions() {
        // a vs b: 2 transitions + 1 transversion over 10 sites -> 0.40235947810852507
        close(k2p(b"AAAACCGGTT", b"AGCACTGGTT"), 0.402_359_478_108_525_07);
    }

    #[test]
    fn k2p_classifies_all_twelve_substitutions() {
        // Each ordered base pair: confirm ts (XOR==2) vs tv classification by
        // checking a one-substitution sequence against the known counts.
        // Transitions: A<->G, C<->T.
        for (x, y) in [(b'A', b'G'), (b'G', b'A'), (b'C', b'T'), (b'T', b'C')] {
            // single transition over 1 site: P=1, Q=0 -> 1-2P-Q = -1 <= 0 -> NaN.
            assert!(k2p(&[x], &[y]).is_nan(), "{}->{} ts", x as char, y as char);
        }
        // Transversions: the 8 purine<->pyrimidine pairs.
        for (x, y) in [
            (b'A', b'C'),
            (b'C', b'A'),
            (b'A', b'T'),
            (b'T', b'A'),
            (b'G', b'C'),
            (b'C', b'G'),
            (b'G', b'T'),
            (b'T', b'G'),
        ] {
            // single transversion over 1 site: P=0, Q=1 -> 1-2Q = -1 <= 0 -> NaN.
            assert!(k2p(&[x], &[y]).is_nan(), "{}->{} tv", x as char, y as char);
        }
    }

    #[test]
    fn k2p_saturation_is_nan() {
        assert!(k2p(b"AAAAAAAA", b"CCCCCCCC").is_nan());
    }

    #[test]
    fn k2p_identical_is_positive_zero() {
        let d = k2p(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(d, 0.0);
        assert!(d.is_sign_positive());
    }

    #[test]
    fn k2p_gap_and_ambiguous_excluded() {
        // N and gaps drop out of the denominator (canonical filter), so a pair
        // differing only at an N column is distance 0 (scikit-bio k2p).
        assert_eq!(k2p(b"ACGTACGTACGT", b"ACGTNCGTACGT"), 0.0);
        // With a gap excluded but a real transition retained, the gap does not
        // dilute the denominator: 1 transition over 7 valid sites (the gap and
        // the last column drop out), matching scikit-bio's pairwise deletion.
        let with_gap = k2p(b"ACGTACG-", b"GCGTACG-");
        let no_gap = k2p(b"ACGTACG", b"GCGTACG");
        close(with_gap, no_gap);
    }

    #[test]
    fn u_equals_t_in_canonical_metrics_but_not_hamming() {
        // K2P maps T and U to the same canonical code, so a T/U column matches.
        assert_eq!(k2p(b"ACGT", b"ACGU"), 0.0);
        // Hamming compares characters literally (like scikit-bio), so T != U.
        close(hamming(b"ACGT", b"ACGU"), 0.25);
    }

    #[test]
    fn metric_from_str() {
        assert_eq!("hamming".parse::<Metric>().unwrap(), Metric::Hamming);
        assert_eq!("jc69".parse::<Metric>().unwrap(), Metric::Jc69);
        assert_eq!("k2p".parse::<Metric>().unwrap(), Metric::K2p);
        assert!("bogus".parse::<Metric>().is_err());
    }
}
