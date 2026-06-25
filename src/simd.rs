//! Stable portable-SIMD canonical-site tallies over pre-encoded code slices
//! (`wide::u8x16`). Codes are 0-3 for A/C/G/T-or-U and 4 for non-canonical
//! (gap/N/ambiguous). A lane is a valid site when neither code is 4.
//!
//! Lane comparisons yield 0xFF/0x00 byte masks; `move_mask` collapses them to a
//! 16-bit word whose `count_ones` popcounts the matched lanes. The scalar tail
//! handles the final `< 16` positions. Results are bit-identical to the scalar
//! `tally_codes`.

use std::ops::{BitAnd, BitOr, BitXor};

use wide::u8x16;

const LANES: usize = 16;
const NON_CANON: u8x16 = u8x16::new([4; LANES]);
const ALL_ONES: u8x16 = u8x16::new([0xFF; LANES]);
const TWO: u8x16 = u8x16::new([2; LANES]);
const ONE: u8x16 = u8x16::new([1; LANES]);

#[inline]
fn load(chunk: &[u8]) -> u8x16 {
    let mut buf = [0u8; LANES];
    buf.copy_from_slice(chunk);
    u8x16::new(buf)
}

/// `wide` 0.7's `u8x16` has no `Not`; complement a byte mask via XOR with all
/// ones (0xFF -> 0x00 and back).
#[inline]
fn not(m: u8x16) -> u8x16 {
    m.bitxor(ALL_ONES)
}

#[inline]
fn popcount(mask: u8x16) -> u64 {
    u64::from((mask.move_mask() as u16).count_ones())
}

/// 0xFF on lanes where both codes are canonical (neither is 4).
#[inline]
fn valid_mask(ca: u8x16, cb: u8x16) -> u8x16 {
    not(ca.cmp_eq(NON_CANON).bitor(cb.cmp_eq(NON_CANON)))
}

/// `(valid, differing)` over canonical sites — the JC69 inputs. `valid` counts
/// positions where both codes are canonical; `differing` counts those that also
/// differ.
#[inline]
#[must_use]
pub fn tally_jc69(a: &[u8], b: &[u8]) -> (u64, u64) {
    debug_assert_eq!(a.len(), b.len());
    let n = a.len();
    let mut valid = 0u64;
    let mut diff = 0u64;

    let mut i = 0;
    while i + LANES <= n {
        let ca = load(&a[i..i + LANES]);
        let cb = load(&b[i..i + LANES]);
        let vmask = valid_mask(ca, cb);
        let dmask = not(ca.cmp_eq(cb)).bitand(vmask);
        valid += popcount(vmask);
        diff += popcount(dmask);
        i += LANES;
    }
    for (&cx, &cy) in a[i..].iter().zip(&b[i..]) {
        let both = cx != 4 && cy != 4;
        valid += u64::from(both);
        diff += u64::from(both && cx != cy);
    }
    (valid, diff)
}

/// `(valid, transitions, transversions)` over canonical sites — the K2P inputs.
/// A transition has code XOR == 2; a transversion has an odd XOR. Both are
/// counted only on valid lanes.
#[inline]
#[must_use]
pub fn tally_k2p(a: &[u8], b: &[u8]) -> (u64, u64, u64) {
    debug_assert_eq!(a.len(), b.len());
    let n = a.len();
    let mut valid = 0u64;
    let mut ts = 0u64;
    let mut tv = 0u64;

    let mut i = 0;
    while i + LANES <= n {
        let ca = load(&a[i..i + LANES]);
        let cb = load(&b[i..i + LANES]);
        let vmask = valid_mask(ca, cb);
        let sub = ca.bitxor(cb);
        let ts_mask = sub.cmp_eq(TWO).bitand(vmask);
        let tv_mask = sub.bitand(ONE).cmp_eq(ONE).bitand(vmask);
        valid += popcount(vmask);
        ts += popcount(ts_mask);
        tv += popcount(tv_mask);
        i += LANES;
    }
    for (&cx, &cy) in a[i..].iter().zip(&b[i..]) {
        if cx == 4 || cy == 4 {
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
    (valid, ts, tv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::{encode, tally_codes};

    fn check(a: &[u8], b: &[u8]) {
        let ea = encode(a);
        let eb = encode(b);
        let (sv, sts, stv) = tally_codes(&ea, &eb);
        let (jv, jd) = tally_jc69(&ea, &eb);
        let (kv, kts, ktv) = tally_k2p(&ea, &eb);
        assert_eq!(jv, sv, "jc69 valid");
        assert_eq!(jd, sts + stv, "jc69 differing");
        assert_eq!((kv, kts, ktv), (sv, sts, stv), "k2p tally");
    }

    #[test]
    fn matches_scalar_across_lengths_and_chars() {
        // Lengths spanning the SIMD body, the tail, and exact multiples of 16.
        check(b"ACGT", b"AGGT");
        check(b"ACGTACGTACGTACGT", b"AGGTACG-ACGTNCGT");
        check(
            b"ACGTACGTACGTACGTACGTNNNN-.gtACGT",
            b"GCATACGTACGTACGTAC-TACGTACGTACGU",
        );
        check(b"", b"");
        check(b"AC-N.tU", b"AGGT.NU");
    }

    #[test]
    fn all_canonical_chars_classify() {
        check(b"AAAAAAAAAAAAAAAAAAAA", b"ACGTACGTACGTACGTACGT");
    }
}
