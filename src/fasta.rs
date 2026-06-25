//! Minimal aligned-FASTA reader. All records must share one alignment length;
//! a ragged alignment is a fail-loud error, not a silent truncation.

use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

pub struct Alignment {
    pub ids: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
}

pub fn read(reader: &mut dyn BufRead) -> Result<Alignment> {
    let mut ids = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut cur: Option<Vec<u8>> = None;

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(seq) = cur.take() {
                seqs.push(seq);
            }
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            if id.is_empty() {
                return Err(RsomicsError::InvalidInput(
                    "FASTA header has no sequence id".into(),
                ));
            }
            ids.push(id);
            cur = Some(Vec::new());
        } else {
            let seq = cur.as_mut().ok_or_else(|| {
                RsomicsError::InvalidInput("sequence data before any '>' header".into())
            })?;
            seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(seq) = cur.take() {
        seqs.push(seq);
    }

    if ids.is_empty() {
        return Err(RsomicsError::InvalidInput("no FASTA records found".into()));
    }
    if ids.len() < 2 {
        return Err(RsomicsError::InvalidInput(format!(
            "need at least 2 sequences for a distance matrix, found {}",
            ids.len()
        )));
    }

    let len = seqs[0].len();
    for (id, seq) in ids.iter().zip(&seqs) {
        if seq.len() != len {
            return Err(RsomicsError::InvalidInput(format!(
                "sequences are not aligned: {:?} has length {}, expected {len}",
                id,
                seq.len()
            )));
        }
    }

    Ok(Alignment { ids, seqs })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parses_aligned_records() {
        let mut c = Cursor::new(">a\nACGT\n>b desc\nAGGT\n");
        let aln = read(&mut c).unwrap();
        assert_eq!(aln.ids, vec!["a", "b"]);
        assert_eq!(aln.seqs[0], b"ACGT");
        assert_eq!(aln.seqs[1], b"AGGT");
    }

    #[test]
    fn wraps_multiline_sequences() {
        let mut c = Cursor::new(">a\nAC\nGT\n>b\nAG\nGT\n");
        let aln = read(&mut c).unwrap();
        assert_eq!(aln.seqs[0], b"ACGT");
    }

    #[test]
    fn ragged_alignment_fails() {
        let mut c = Cursor::new(">a\nACGT\n>b\nACG\n");
        assert!(read(&mut c).is_err());
    }

    #[test]
    fn single_record_fails() {
        let mut c = Cursor::new(">a\nACGT\n");
        assert!(read(&mut c).is_err());
    }

    #[test]
    fn lowercase_is_uppercased() {
        let mut c = Cursor::new(">a\nacgt\n>b\nacgt\n");
        let aln = read(&mut c).unwrap();
        assert_eq!(aln.seqs[0], b"ACGT");
    }
}
