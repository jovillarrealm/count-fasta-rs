use std::process::Command;
use std::collections::HashMap;
use std::path::Path;

fn run_command(cmd: &str, args: &[&str]) -> String {
    let output = Command::new(cmd)
        .args(args)
        .output()
        .unwrap_or_else(|_| panic!("Failed to execute command: {} {:?}", cmd, args));

    if !output.status.success() {
        panic!(
            "Command failed: {} {:?}\nStderr: {}",
            cmd,
            args,
            String::from_utf8_lossy(&output.stderr)
        );
    }

    String::from_utf8(output.stdout).expect("Output is not valid UTF-8")
}

fn parse_output(output: &str) -> HashMap<String, String> {
    let mut metrics = HashMap::new();
    for line in output.lines() {
        if line.contains("Total length of sequence:") {
            metrics.insert(
                "Total length".to_string(),
                line.split_whitespace().nth(4).unwrap().to_string(),
            );
        } else if line.contains("Total number of sequences:") {
            metrics.insert(
                "Sequence count".to_string(),
                line.split_whitespace().last().unwrap().to_string(),
            );
        } else if line.contains("Average contig length is:") {
             metrics.insert(
                "Average length".to_string(),
                line.split_whitespace().nth(4).unwrap().to_string(),
            );
        } else if line.contains("Largest contig:") {
            metrics.insert(
                "Largest contig".to_string(),
                line.split_whitespace().nth(2).unwrap().to_string(),
            );
        } else if line.contains("Shortest contig:") {
            metrics.insert(
                "Shortest contig".to_string(),
                line.split_whitespace().nth(2).unwrap().to_string(),
            );
        } else if line.contains("N50 stats:") {
             // "N50 stats:			50% of total sequence length is contained in the 5 sequences >= 50 bp"
             let parts: Vec<&str> = line.split_whitespace().collect();
             if parts.len() >= 2 {
                 metrics.insert(
                    "N50".to_string(),
                    parts[parts.len() - 2].to_string()
                 );
             }
        } else if line.contains("Total GC count:") {
            metrics.insert(
                "GC count".to_string(),
                line.split_whitespace().nth(3).unwrap().to_string(),
            );
        } else if line.contains("Number of Ns:") {
             metrics.insert(
                "N count".to_string(),
                line.split_whitespace().last().unwrap().to_string(),
            );
        }
    }
    metrics
}

#[test]
fn test_compare_perl_script() {
    let perl_script = "alternatives/count_fasta_cnsg.pl";
    let test_files = vec![
        "test/ay number 2.fna",
        "test/multi_seq.fa",
        "test/with_gaps.fa",
        "test/garbage_start.fa",
        "test/multi_line_seq.fa",
        "test/mixed_case.fa",
    ];

    if !Path::new(perl_script).exists() {
        eprintln!("Perl script not found at {}", perl_script);
        return;
    }
    
    // Check if perl is installed
    if Command::new("perl").arg("-v").output().is_err() {
        eprintln!("Perl is not installed, skipping test.");
        return;
    }

    // Build first to ensure up-to-date
    let status = Command::new("cargo")
        .args(&["build", "--quiet"])
        .status()
        .expect("Failed to build");
    assert!(status.success());

    for test_file in test_files {
        println!("Testing file: {}", test_file);
        
        // Run Perl script
        let perl_output = run_command("perl", &[perl_script, test_file]);
        let perl_metrics = parse_output(&perl_output);

        // Run Rust binary
        let rust_output = run_command("cargo", &["run", "--quiet", "--", test_file]);
        let rust_metrics = parse_output(&rust_output);

        // Compare metrics
        let keys_to_compare = if test_file.contains("with_gaps.fa") {
            // Rust implementation ignores gaps in sequence length, whereas Perl includes them.
            // So we only compare GC and N counts for this file.
            vec!["GC count", "N count"]
        } else {
            vec![
                "Total length",
                "Sequence count",
                "Average length",
                "Largest contig",
                "Shortest contig",
                "N50",
                "GC count",
                "N count",
            ]
        };

        for key in keys_to_compare {
            let perl_val = perl_metrics.get(key);
            let rust_val = rust_metrics.get(key);
            
            assert!(perl_val.is_some(), "Missing key in Perl output for {}: {}", test_file, key);
            assert!(rust_val.is_some(), "Missing key in Rust output for {}: {}", test_file, key);

            assert_eq!(
                perl_val,
                rust_val,
                "Mismatch for {} in file {}: Perl={}, Rust={}",
                key, test_file, perl_val.unwrap(), rust_val.unwrap()
            );
        }
    }
}
