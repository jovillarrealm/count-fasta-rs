// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

//! # count-fasta-rs
//!
//! A high-performance stream processor for FASTA files (including compressed formats) 
//! designed with the following architecture:
//! - **Low memory usage**: Processes files in chunks to minimize footprint.
//! - **Zero-copy reads**: Leverages memory mapping and buffered I/O to avoid unnecessary data duplication.
//! - **One file per core**: Utilizes parallel processing with Rayon, scaling efficiently across available CPUs.
//! - **Efficient I/O**: Optimizes OS-level read-ahead and sequential access patterns.
//! - **SIMD optimizations**: Employs AVX2 instructions for rapid sequence analysis and statistics calculation.

use clap::Parser;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::env;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

mod process_files;
mod simd;

extern crate num_cpus;

fn determine_buffer_size() -> usize {
    const DEFAULT_SIZE: usize = 2 * 1024 * 1024; // Default to 1MB
    const MAX_SIZE: usize = 5 * 1024 * 1024; // Cap at 10MB
    match env::var("BUFFER_SIZE") {
        Ok(val) => val.parse().unwrap_or(DEFAULT_SIZE).min(MAX_SIZE),
        Err(_) => DEFAULT_SIZE,
    }
}

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about = "A high-performance stream processor for FASTA files.",
    long_about = "🚀 A high-performance, SIMD-accelerated stream processor for FASTA files.

This tool calculates assembly statistics (N50, GC%, total length, etc.) for FASTA files, supporting a wide range of compression formats. It leverages multi-threading (Rayon) and vectorization (AVX2/AVX512/NEON) to process data at maximum speed.

SUPPORTED FORMATS:
  • Uncompressed: .fa, .fasta, .fna
  • Compressed:   .gz, .bgz, .bgzip (Block GZIP), .xz, .bz2, .naf (Nucleotide Archive)
  • Archives:     .zip (processes all valid FASTA files inside)

TUTORIAL & EXAMPLES:

1. Basic Usage (Standard Output)
   Calculate stats for a single file and print to stdout:
     $ count-fasta-rs genome.fna

2. Processing Multiple Files
   Use wildcards to process multiple files at once:
     $ count-fasta-rs *.fa.gz

3. Directory Processing
   Process all valid files within a directory (non-recursive). You can specify multiple directories:
     $ count-fasta-rs -d ./genomes -d ./more_genomes

4. Saving Results to CSV
   Append results to a CSV file for easy analysis in Excel/Pandas:
     $ count-fasta-rs -c results.csv -d ./genomes
   (Note: If 'results.csv' exists, new rows are appended. If not, it's created with a header.)

5. Performance Tuning
   - Threads: By default, it uses all available cores. Limit this with -t:
     $ count-fasta-rs -t 4 genome.fna
   - SIMD: If you encounter issues or want to compare scalar performance:
     $ count-fasta-rs --no-simd genome.fna

NOTES:
  - Gzip (.gz) files are assumed to be standard gzip. For random access optimized BGZIP, use .bgz or .bgzip extensions if possible, though .gz works for sequential reading.
  - The tool uses zero-copy reading where possible (mmap) to keep memory usage low, even for huge files."
)]
struct Args {
    /// Path to the CSV file to be created or appended to.
    ///
    /// If the file does not exist, it will be created with a header.
    /// If it exists, new results will be appended.
    #[clap(short, long, value_hint = clap::ValueHint::FilePath)]
    csv: Option<String>,

    /// Directory to process (non-recursive).
    ///
    /// The program will process all valid FASTA files found in these directories.
    #[clap(short, long, value_hint = clap::ValueHint::DirPath)]
    directory: Vec<String>,

    /// Number of threads to use.
    ///
    /// If not specified, the program will automatically determine the number of threads based on
    /// available CPUs and file count.
    #[clap(short, long)]
    threads: Option<usize>,

    /// Use legacy output format.
    ///
    /// Useful for debugging and testing purposes to match older output styles.
    #[clap(short, long)]
    legacy: bool,

    /// FASTA file(s) to process.
    ///
    /// Supports wildcards. Inside a zip file, only .fa, .fasta, and .fna files will be processed.
    /// Gzip (.gz) files are assumed to be standard gzip; bgzip files should ideally use .bgz or .bgzip.
    #[clap(name = "FASTA FILE", value_hint = clap::ValueHint::FilePath)]
    files: Vec<String>,

    /// Disable SIMD optimizations (force scalar fallback).
    ///
    /// Useful for debugging or if SIMD causes issues on specific hardware.
    #[clap(long)]
    no_simd: bool,
}

fn main() {
    let args = Args::parse();
    let mut files_to_process = Vec::new();

    for dir in args.directory {
        match get_fasta_files_from_directory(&dir) {
            Ok(files) => files_to_process.extend(files),
            Err(e) => {
                eprintln!("Error reading directory '{}': {}", dir, e);
                std::process::exit(1);
            }
        }
    }
    files_to_process.extend(args.files.into_iter().map(PathBuf::from));

    let results = process_files(files_to_process, args.threads, args.no_simd);

    if let Some(csv_file) = args.csv {
        append_to_csv(&results, &csv_file).expect("Failed to write CSV");
    } else {
        for result in results {
            print_results(&result, args.legacy);
        }
    }
}

fn get_fasta_files_from_directory(dir: &str) -> std::io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file() && path.extension().and_then(|s| s.to_str()).is_some_and(|ext| {
            process_files::VALID_FILES.contains(&ext) || process_files::VALID_COMPRESSION.contains(&ext)
        }) {
            files.push(path);
        }
    }
    Ok(files)
}

fn process_files(
    files: Vec<PathBuf>,
    threads: Option<usize>,
    no_simd: bool,
) -> Vec<process_files::AnalysisResults> {
    let buffer_size = determine_buffer_size();
    let available_threads = determine_threads(&files, threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(available_threads)
        .build()
        .unwrap();
    pool.install(|| {
        files
            .par_iter()
            .flat_map(|file| {
                let ext = file.extension().and_then(|e| e.to_str()).unwrap_or("");
                let res = match ext {
                    "gz" => process_files::process_gz_file(file, buffer_size, no_simd),
                    "zip" => process_files::process_zip_file(file, buffer_size, no_simd),
                    "xz" => process_files::process_xz_file(file, buffer_size, no_simd),
                    "bz2" => process_files::process_bz2_file(file, buffer_size, no_simd),
                    "bgz" | "bgzip" => {
                        process_files::process_bgzip_file(file, buffer_size, no_simd)
                    }
                    "naf" => process_files::process_naf_file(file, no_simd),
                    _ if process_files::VALID_FILES.contains(&ext) => {
                        process_files::process_fasta_file(file, buffer_size, no_simd)
                    }
                    _ => Ok(Vec::new()),
                };
                match res {
                    Ok(v) => v,
                    Err(e) => {
                        eprintln!("Error processing file {:?}: {}", file, e);
                        Vec::new()
                    }
                }
            })
            .collect()
    })
}

fn determine_threads(files: &[PathBuf], threads: Option<usize>) -> usize {
    let available_threads;
    if let Some(threads) = threads {
        available_threads = threads;
    } else {
        let usable_threads_logical = (num_cpus::get() as f32 * 0.9).round() as usize;
        let usable_physical_threads = (num_cpus::get_physical() as f32 * 0.75).round() as usize;
        let usable_threads = max(usable_threads_logical, usable_physical_threads);
        available_threads = min(usable_threads, files.len());
    }
    available_threads
}

fn print_results(results: &process_files::AnalysisResults, legacy: bool) {
    if !legacy {
        println!("\nFile name:\t{} ", results.filename);
    } else {
        println!();
    }
    println!("Total length of sequence:\t{} bp", results.total_length);
    println!("Total number of sequences:\t{}", results.sequence_count);
    let avg_len = if results.sequence_count > 0 {
        results.total_length / results.sequence_count
    } else {
        0
    };
    println!(
        "Average contig length is:\t{} bp",
        avg_len
    );
    println!("Largest contig:\t\t{} bp", results.largest_contig);
    println!("Shortest contig:\t\t{} bp", results.shortest_contig);
    println!(
        "N25 stats:\t\t\t25% of total sequence length is contained in the {} sequences >= {} bp",
        results.n25_sequence_count, results.n25
    );
    println!(
        "N50 stats:\t\t\t50% of total sequence length is contained in the {} sequences >= {} bp",
        results.n50_sequence_count, results.n50
    );
    println!(
        "N75 stats:\t\t\t75% of total sequence length is contained in the {} sequences >= {} bp",
        results.n75_sequence_count, results.n75
    );
    println!("Total GC count:\t\t\t{} bp", results.gc_count);
    println!(
        "GC %:\t\t\t\t{:.2} %",
        (results.gc_count as f64 / results.total_length as f64) * 100.0
    );
    println!("Number of Ns:\t\t\t{}", results.n_count);
    println!(
        "Ns %:\t\t\t\t{:.2} %",
        (results.n_count as f64 / results.total_length as f64) * 100.0
    );
}

fn append_to_csv(results: &[process_files::AnalysisResults], csv_filename: &str) -> io::Result<()> {
    let csv_exists = Path::new(csv_filename).exists();
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_filename)?;
    let mut writer = std::io::BufWriter::new(file);

    if !csv_exists {
        let header = "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n";
        writer.write_all(header.as_bytes())?;
    }

    for result in results {
        let avg_len = if result.sequence_count > 0 {
            (result.total_length as f64 / result.sequence_count as f64).round() as usize
        } else {
            0
        };
        let gc_pct = if result.total_length > 0 {
            (result.gc_count as f64 / result.total_length as f64) * 100.0
        } else {
            0.0
        };
        let n_pct = if result.total_length > 0 {
            (result.n_count as f64 / result.total_length as f64) * 100.0
        } else {
            0.0
        };

        writeln!(
            writer,
            "{};{};{};{};{};{};{};{:.7};{};{:.7}",
            result.filename,
            result.total_length,
            result.sequence_count,
            avg_len,
            result.largest_contig,
            result.shortest_contig,
            result.n50,
            gc_pct,
            result.n_count,
            n_pct,
        )?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use std::iter::zip as zip_things;
    #[test]
    fn it_works() {
        let mut files_to_process = Vec::new();

        if let Ok(files) = get_fasta_files_from_directory(&"./test/") {
            files_to_process.extend(files);
        }

        let results = process_files(files_to_process, None, false);

        let csv_file = "test/attempt.csv";
        if Path::new(csv_file).exists() {
            let _ = fs::remove_file(csv_file);
        }

        append_to_csv(&results, &csv_file).expect("Failed to write CSV");
        let mut thing: Vec<String> = fs::read_to_string("test/test.csv")
            .unwrap()
            .lines()
            .map(String::from) // make each slice into a string
            .collect();
        let mut compare: Vec<String> = fs::read_to_string(csv_file)
            .unwrap()
            .lines()
            .map(String::from) // make each slice into a string
            .collect();
        thing.sort();
        compare.sort();
        for (thing_i, compare_i) in zip_things(thing, compare) {
            assert_eq!(thing_i, compare_i);
        }
        let _ = fs::remove_file(csv_file);
    }
}
