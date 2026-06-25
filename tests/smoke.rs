use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-seq-dist"))
}

fn fixture() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/aligned.fasta")
}

#[test]
fn default_metric_is_hamming() {
    let with = Command::new(ours())
        .arg(fixture())
        .args(["--metric", "hamming"])
        .output()
        .unwrap();
    let without = Command::new(ours()).arg(fixture()).output().unwrap();
    assert_eq!(with.stdout, without.stdout);
}

#[test]
fn output_is_square_with_id_header() {
    let out = Command::new(ours())
        .arg(fixture())
        .args(["--metric", "jc69"])
        .output()
        .unwrap();
    let text = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let n = lines[0].split('\t').count();
    assert_eq!(lines.len(), n + 1, "header + n rows");
    for row in &lines[1..] {
        assert_eq!(row.split('\t').count(), n, "each row has n columns");
    }
}

#[test]
fn unknown_metric_rejected() {
    let out = Command::new(ours())
        .arg(fixture())
        .args(["--metric", "bogus"])
        .output()
        .unwrap();
    assert!(!out.status.success());
}

#[test]
fn output_to_file() {
    let dir = tempdir();
    let path = dir.join("out.tsv");
    let out = Command::new(ours())
        .arg(fixture())
        .args(["--metric", "k2p", "-o"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(out.status.success());
    let written = std::fs::read_to_string(&path).unwrap();
    assert!(written.starts_with("s1\ts2"));
}

fn tempdir() -> std::path::PathBuf {
    let base =
        std::env::var("TMPDIR").unwrap_or_else(|_| "/Volumes/KIOXIA/Developments/tmp".to_string());
    let d = std::path::PathBuf::from(base).join(format!("seq-dist-smoke-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}
