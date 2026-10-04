use std::process::Command;

#[test]
fn test_cli_min_length_all_excluded() {
    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["-m", "10", "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Total length of sequence:\t0 bp"));
    assert!(stdout.contains("Total number of sequences:\t0"));
}

#[test]
fn test_cli_min_length_partial() {
    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["-m", "7", "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Total length of sequence:\t16 bp"));
    assert!(stdout.contains("Total number of sequences:\t2"));
    assert!(stdout.contains("Largest contig:\t\t8 bp"));
    assert!(stdout.contains("Shortest contig:\t\t8 bp"));
}
