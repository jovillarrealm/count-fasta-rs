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

#[test]
fn test_cli_csv_rfc4180() {
    let temp_csv = "test/test_rfc4180.csv";
    if std::path::Path::new(temp_csv).exists() {
        let _ = std::fs::remove_file(temp_csv);
    }

    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["-c", temp_csv, "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let content = std::fs::read_to_string(temp_csv).expect("Failed to read CSV");
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "filename,total_length,sequence_count,average_length,largest_contig,shortest_contig,n25,l25,n50,l50,n75,l75,n90,l90,gc_count,gc_pct,n_count,n_pct");
    let _ = std::fs::remove_file(temp_csv);
}

#[test]
fn test_cli_tsv_export() {
    let temp_tsv = "test/test_export.tsv";
    if std::path::Path::new(temp_tsv).exists() {
        let _ = std::fs::remove_file(temp_tsv);
    }

    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["--tsv", temp_tsv, "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let content = std::fs::read_to_string(temp_tsv).expect("Failed to read TSV");
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "filename\ttotal_length\tsequence_count\taverage_length\tlargest_contig\tshortest_contig\tn25\tl25\tn50\tl50\tn75\tl75\tn90\tl90\tgc_count\tgc_pct\tn_count\tn_pct");
    let _ = std::fs::remove_file(temp_tsv);
}

#[test]
fn test_cli_json_export() {
    let temp_json = "test/test_export.json";
    if std::path::Path::new(temp_json).exists() {
        let _ = std::fs::remove_file(temp_json);
    }

    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["--json", temp_json, "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let content = std::fs::read_to_string(temp_json).expect("Failed to read JSON");
    let v: serde_json::Value = serde_json::from_str(&content).expect("Valid JSON");
    assert!(v.is_array());
    assert_eq!(v[0]["total_length"], 24);
    assert_eq!(v[0]["sequence_count"], 4);
    let _ = std::fs::remove_file(temp_json);
}

#[test]
fn test_cli_legacy_csv_export() {
    let temp_csv = "test/test_legacy.csv";
    if std::path::Path::new(temp_csv).exists() {
        let _ = std::fs::remove_file(temp_csv);
    }

    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["-c", temp_csv, "--legacy", "test/multi_seq.fa"])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let content = std::fs::read_to_string(temp_csv).expect("Failed to read legacy CSV");
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines[0],
        "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage"
    );
    assert!(lines[1].contains("multi_seq.fa;24;4;6;8;4;8;"));
    let _ = std::fs::remove_file(temp_csv);
}

#[test]
fn test_cli_zip_archive_provenance() {
    use std::io::Write;
    let zip_path = "test/test_archive.zip";
    let temp_json = "test/test_archive.json";

    let file = std::fs::File::create(zip_path).expect("Failed to create zip");
    let mut zip = zip::ZipWriter::new(file);
    let options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Stored);
    zip.start_file("nested/test_seq.fa", options).unwrap();
    zip.write_all(b">seq1\nATGCATGC\n").unwrap();
    zip.finish().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_count-fasta-rs"))
        .args(&["--json", temp_json, zip_path])
        .output()
        .expect("Failed to execute binary");

    assert!(output.status.success());
    let content = std::fs::read_to_string(temp_json).expect("Failed to read JSON");
    let v: serde_json::Value = serde_json::from_str(&content).expect("Valid JSON");
    assert_eq!(v[0]["filename"], "test_archive.zip#nested/test_seq.fa");
    assert_eq!(v[0]["total_length"], 8);
    assert_eq!(v[0]["sequence_count"], 1);

    let _ = std::fs::remove_file(zip_path);
    let _ = std::fs::remove_file(temp_json);
}

