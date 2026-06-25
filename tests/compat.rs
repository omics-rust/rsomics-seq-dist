//! Always-run compat: the binary's distance-matrix TSV must reproduce the
//! committed scikit-bio 0.7.2 goldens byte-for-byte. No scikit-bio at test time —
//! the goldens are checked in under tests/golden/.

use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-seq-dist"))
}

fn golden(name: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("tests/golden/{name}"))
}

fn run_metric(metric: &str) -> String {
    let out = Command::new(ours())
        .arg(golden("aligned.fasta"))
        .args(["--metric", metric])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn check(metric: &str) {
    let got = run_metric(metric);
    let want = std::fs::read_to_string(golden(&format!("{metric}.golden.tsv"))).unwrap();
    assert_eq!(got, want, "{metric} output diverged from scikit-bio golden");
}

#[test]
fn hamming_matches_golden() {
    check("hamming");
}

#[test]
fn jc69_matches_golden() {
    check("jc69");
}

#[test]
fn k2p_matches_golden() {
    check("k2p");
}

#[test]
fn stdin_matches_file() {
    use std::io::Write;
    let fasta = std::fs::read(golden("aligned.fasta")).unwrap();
    let mut child = Command::new(ours())
        .args(["-", "--metric", "k2p"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&fasta).unwrap();
    let out = child.wait_with_output().unwrap();
    let got = String::from_utf8(out.stdout).unwrap();
    let want = std::fs::read_to_string(golden("k2p.golden.tsv")).unwrap();
    assert_eq!(got, want);
}
