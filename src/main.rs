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

use clap::{CommandFactory, FromArgMatches, Parser};
use rayon::prelude::*;
use std::cmp::min;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

mod process_files;
mod simd;

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

    /// Path to the TSV file to be created or appended to.
    ///
    /// If the file does not exist, it will be created with a header.
    /// If it exists, new results will be appended.
    #[clap(long, value_hint = clap::ValueHint::FilePath)]
    tsv: Option<String>,

    /// Path to the JSON file to write results to.
    #[clap(long, value_hint = clap::ValueHint::FilePath)]
    json: Option<String>,

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

    /// Minimum contig length threshold in bp.
    ///
    /// Contigs shorter than this length are excluded from assembly metrics.
    #[clap(short = 'm', long = "min-length", default_value = "0")]
    min_length: usize,

    /// Disable SIMD optimizations (force scalar fallback).
    ///
    /// Useful for debugging or if SIMD causes issues on specific hardware.
    #[clap(long)]
    no_simd: bool,
}

fn main() {
    let cmd = Args::command().after_help(format!(
        "SIMD Support: Enabled (std::simd)\nDefault max threads: {}",
        determine_threads(None, None)
    ));

    let matches = cmd.get_matches();
    let args = Args::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    let mut files_to_process = Vec::new();
    let mut directory_error = false;
    for dir in args.directory {
        match get_fasta_files_from_directory(&dir) {
            Ok(files) => files_to_process.extend(files),
            Err(e) => {
                eprintln!("Error reading directory '{}': {}", dir, e);
                directory_error = true;
            }
        }
    }
    files_to_process.extend(args.files.into_iter().map(PathBuf::from));

    if files_to_process.is_empty() && !directory_error {
        files_to_process.push(PathBuf::from("-"));
    }

    let (results, processing_error) =
        process_files(files_to_process, args.threads, args.no_simd, args.min_length);
    let mut exit_error = directory_error || processing_error;

    let has_file_output = args.csv.is_some() || args.tsv.is_some() || args.json.is_some();

    if let Some(ref csv_file) = args.csv {
        if let Err(e) = append_to_csv(&results, csv_file, args.legacy) {
            eprintln!("Failed to write CSV file '{}': {}", csv_file, e);
            exit_error = true;
        }
    }

    if let Some(ref tsv_file) = args.tsv {
        if let Err(e) = append_to_tsv(&results, tsv_file) {
            eprintln!("Failed to write TSV file '{}': {}", tsv_file, e);
            exit_error = true;
        }
    }

    if let Some(ref json_file) = args.json {
        if let Err(e) = write_to_json(&results, json_file) {
            eprintln!("Failed to write JSON file '{}': {}", json_file, e);
            exit_error = true;
        }
    }

    if !has_file_output {
        for result in &results {
            print_results(result, args.legacy);
        }
    }

    if exit_error {
        std::process::exit(1);
    }
}

fn get_fasta_files_from_directory(dir: &str) -> std::io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file()
            && process_files::FileFormat::from_path(&path) != process_files::FileFormat::Unknown
        {
            files.push(path);
        }
    }
    Ok(files)
}

fn process_files(
    files: Vec<PathBuf>,
    threads: Option<usize>,
    no_simd: bool,
    min_length: usize,
) -> (Vec<process_files::AnalysisResults>, bool) {
    let available_threads = determine_threads(Some(&files), threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(available_threads)
        .build()
        .unwrap();
    let had_error = std::sync::atomic::AtomicBool::new(false);
    let results: Vec<process_files::AnalysisResults> = pool.install(|| {
        files
            .par_iter()
            .flat_map(
                |file| match process_files::process_any_file(file, no_simd, min_length) {
                    Ok(v) => v,
                    Err(e) => {
                        eprintln!("Error processing file {:?}: {}", file, e);
                        had_error.store(true, std::sync::atomic::Ordering::Relaxed);
                        Vec::new()
                    }
                },
            )
            .collect()
    });
    (results, had_error.load(std::sync::atomic::Ordering::Relaxed))
}

fn determine_threads(files: Option<&[PathBuf]>, threads: Option<usize>) -> usize {
    let available_threads;
    if let Some(threads) = threads {
        available_threads = threads;
    } else {
        let usable_threads_logical = (num_cpus::get() as f32 * 0.9).round() as usize;
        let usable_physical_threads = (num_cpus::get_physical() as f32 * 0.75).round() as usize;
        let usable_threads = min(usable_threads_logical, usable_physical_threads);

        available_threads = min(usable_threads, files.map_or(usize::MAX, |f| f.len()));
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
    let avg_len = results
        .total_length
        .checked_div(results.sequence_count)
        .unwrap_or(0);
    println!("Average contig length is:\t{} bp", avg_len);
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
    println!(
        "N90 stats:\t\t\t90% of total sequence length is contained in the {} sequences >= {} bp",
        results.n90_sequence_count, results.n90
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

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AssemblyReport {
    pub filename: String,
    pub total_length: usize,
    pub sequence_count: usize,
    pub average_length: usize,
    pub largest_contig: usize,
    pub shortest_contig: usize,
    pub n25: usize,
    pub l25: usize,
    pub n50: usize,
    pub l50: usize,
    pub n75: usize,
    pub l75: usize,
    pub n90: usize,
    pub l90: usize,
    pub gc_count: usize,
    pub gc_pct: f64,
    pub n_count: usize,
    pub n_pct: f64,
}

impl From<&process_files::AnalysisResults> for AssemblyReport {
    fn from(result: &process_files::AnalysisResults) -> Self {
        let average_length = if result.sequence_count > 0 {
            (result.total_length as f64 / result.sequence_count as f64).round() as usize
        } else {
            0
        };
        let gc_pct = if result.total_length > 0 {
            ((result.gc_count as f64 / result.total_length as f64) * 10000.0).round() / 100.0
        } else {
            0.0
        };
        let n_pct = if result.total_length > 0 {
            ((result.n_count as f64 / result.total_length as f64) * 10000.0).round() / 100.0
        } else {
            0.0
        };

        Self {
            filename: result.filename.clone(),
            total_length: result.total_length,
            sequence_count: result.sequence_count,
            average_length,
            largest_contig: result.largest_contig,
            shortest_contig: result.shortest_contig,
            n25: result.n25,
            l25: result.l25,
            n50: result.n50,
            l50: result.l50,
            n75: result.n75,
            l75: result.l75,
            n90: result.n90,
            l90: result.l90,
            gc_count: result.gc_count,
            gc_pct,
            n_count: result.n_count,
            n_pct,
        }
    }
}

fn escape_csv(val: &str) -> String {
    if val.contains(',') || val.contains('"') || val.contains('\n') || val.contains('\r') {
        format!("\"{}\"", val.replace('"', "\"\""))
    } else {
        val.to_string()
    }
}

fn append_to_csv(
    results: &[process_files::AnalysisResults],
    csv_filename: &str,
    legacy: bool,
) -> io::Result<()> {
    let csv_exists = Path::new(csv_filename).exists();
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_filename)?;
    let mut writer = std::io::BufWriter::new(file);

    if !csv_exists {
        if legacy {
            let header = "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n";
            writer.write_all(header.as_bytes())?;
        } else {
            let header = "filename,total_length,sequence_count,average_length,largest_contig,shortest_contig,n25,l25,n50,l50,n75,l75,n90,l90,gc_count,gc_pct,n_count,n_pct\n";
            writer.write_all(header.as_bytes())?;
        }
    }

    if legacy {
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
    } else {
        for result in results {
            let report = AssemblyReport::from(result);
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.2},{},{:.2}",
                escape_csv(&report.filename),
                report.total_length,
                report.sequence_count,
                report.average_length,
                report.largest_contig,
                report.shortest_contig,
                report.n25,
                report.l25,
                report.n50,
                report.l50,
                report.n75,
                report.l75,
                report.n90,
                report.l90,
                report.gc_count,
                report.gc_pct,
                report.n_count,
                report.n_pct,
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn append_to_tsv(
    results: &[process_files::AnalysisResults],
    tsv_filename: &str,
) -> io::Result<()> {
    let tsv_exists = Path::new(tsv_filename).exists();
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(tsv_filename)?;
    let mut writer = std::io::BufWriter::new(file);

    if !tsv_exists {
        let header = "filename\ttotal_length\tsequence_count\taverage_length\tlargest_contig\tshortest_contig\tn25\tl25\tn50\tl50\tn75\tl75\tn90\tl90\tgc_count\tgc_pct\tn_count\tn_pct\n";
        writer.write_all(header.as_bytes())?;
    }

    for result in results {
        let report = AssemblyReport::from(result);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.2}",
            report.filename,
            report.total_length,
            report.sequence_count,
            report.average_length,
            report.largest_contig,
            report.shortest_contig,
            report.n25,
            report.l25,
            report.n50,
            report.l50,
            report.n75,
            report.l75,
            report.n90,
            report.l90,
            report.gc_count,
            report.gc_pct,
            report.n_count,
            report.n_pct,
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_to_json(
    results: &[process_files::AnalysisResults],
    json_filename: &str,
) -> io::Result<()> {
    let mut reports: Vec<AssemblyReport> = if Path::new(json_filename).exists() {
        if let Ok(file) = std::fs::File::open(json_filename) {
            serde_json::from_reader(file).unwrap_or_default()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    reports.extend(results.iter().map(AssemblyReport::from));

    let file = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(json_filename)?;
    let mut writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, &reports)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
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

        let (results, _) = process_files(files_to_process, None, false, 0);

        let csv_file = "test/attempt.csv";
        if Path::new(csv_file).exists() {
            let _ = fs::remove_file(csv_file);
        }

        append_to_csv(&results, &csv_file, true).expect("Failed to write CSV");
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
