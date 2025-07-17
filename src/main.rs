// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

use clap::Parser;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::env;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

mod process_files;

extern crate bytecount;
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
    long_about = "Calculates a dir or a FASTA_FILE and prints its output to stdout.
Or to a csv file. If the csv file already exists, it appends.

FASTA_FILE can be .fa .fasta .fna .zip .gz .xz .bz2 .bgz .bgzip
    
For example:
    count-fasta-rs -c stats.csv -d path/GENOMIC
    count-fasta-rs -c stats.csv genomic*files*.f*a"
)]
struct Args {
    /// Path to csv to be created. It will not create a directories if they don't exist.
    ///
    /// It will append to the csv file if it already exists.
    #[clap(short, long)]
    csv: Option<String>,

    /// Directory to be processed. Non-recursively.
    ///
    /// The program will process all any FASTA_FILE in the path.
    #[clap(short, long)]
    directory: Option<String>,

    /// Numbers of threads to be used, otherwise the program will decide on its own.
    ///
    /// It will decide based on the number of available logical threads, physical cpus, checking cgroups, and number of files
    /// to be processed. On older machines it probably would default to 1 so it's better to set it manually when running large amounts of data against this.
    #[clap(short, long)]
    threads: Option<usize>,

    /// Legacy output
    ///
    /// For debugging and testing purposes.
    #[clap(short, long)]
    legacy: bool,

    /// FASTA FILE[s] to be processed [wildcards would work here].
    ///
    /// Inside a zip file, only .fa .fasta .fna files will be processed.
    /// gz files are assumed to be compressed using gzip, wrong results can come out of a gz file compressed using bzip
    #[clap(name = "FASTA FILE")]
    files: Vec<String>,
}

fn main() {
    let args = Args::parse();
    let mut files_to_process = Vec::new();

    if let Some(dir) = args.directory {
        if let Ok(files) = get_fasta_files_from_directory(&dir) {
            files_to_process.extend(files);
        }
    }
    files_to_process.extend(args.files.into_iter().map(PathBuf::from));

    let results = process_files(files_to_process, args.threads);

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
        if path.is_file() {
            if let Some(ext) = path.extension().and_then(|ext| ext.to_str()) {
                if process_files::VALID_FILES.contains(&ext) {
                    files.push(path);
                } else if process_files::VALID_COMPRESSION.contains(&ext) {
                    files.push(path);
                }
            }
        }
    }
    Ok(files)
}

fn process_files(
    files: Vec<PathBuf>,
    threads: Option<usize>,
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
                match ext {
                    "gz" => process_files::process_gz_file(file, buffer_size).unwrap_or_default(),
                    "zip" => process_files::process_zip_file(file, buffer_size).unwrap_or_default(),
                    "xz" => process_files::process_xz_file(file, buffer_size).unwrap_or_default(),
                    "bz2" => process_files::process_bz2_file(file, buffer_size).unwrap_or_default(),
                    "bgz" | "bgzip" => {
                        process_files::process_bgzip_file(file, buffer_size).unwrap_or_default()
                    }
                    "naf" => process_files::process_naf_file(file).unwrap_or_default(),
                    _ if process_files::VALID_FILES.contains(&ext) => {
                        process_files::process_fasta_file(file, buffer_size).unwrap_or_default()
                    }
                    _ => Vec::new(),
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
    println!(
        "Average contig length is:\t{} bp",
        results.total_length / results.sequence_count
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
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_filename)?;

    if !csv_exists {
        let header = "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n";
        file.write_all(header.as_bytes())?;
    }

    let mut buffer = String::new();
    for result in results {
        let line = format!(
            "{};{};{};{};{};{};{};{:.7};{};{:.7}\n",
            result.filename,
            result.total_length,
            result.sequence_count,
            (result.total_length as f64 / result.sequence_count as f64).round() as usize,
            result.largest_contig,
            result.shortest_contig,
            result.n50,
            (result.gc_count as f64 / result.total_length as f64) * 100.0,
            result.n_count,
            (result.n_count as f64 / result.total_length as f64) * 100.0,
        );
        buffer.push_str(&line);

        // Write in chunks to avoid holding too much in memory
        if buffer.len() > 128 * 1024 {
            file.write_all(buffer.as_bytes())?;
            buffer.clear();
        }
    }

    if !buffer.is_empty() {
        file.write_all(buffer.as_bytes())?;
    }

    file.flush()?;
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

        let results = process_files(files_to_process, None);

        let csv_file = "test/attempt.csv";

        append_to_csv(&results, &csv_file).expect("Failed to write CSV");
        let mut thing:Vec<String> = fs::read_to_string("test/test.csv")
            .unwrap()
            .lines()
            .map(String::from) // make each slice into a string
            .collect();
        let mut compare:Vec<String> = fs::read_to_string(csv_file)
            .unwrap()
            .lines()
            .map(String::from) // make each slice into a string
            .collect();
        thing.sort();
        compare.sort();
        for (thing_i, compare_i) in zip_things(thing, compare){
            assert_eq!(thing_i, compare_i);
        } 
        let _ = fs::remove_file(csv_file);
    }
}
