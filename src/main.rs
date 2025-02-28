// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

use bzip2::read::BzDecoder;
use clap::Parser;
use flate2::read::GzDecoder;
use noodles_bgzf as bgzf;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use xz2::read::XzDecoder;
use zip::read::ZipArchive;
extern crate bytecount;
extern crate num_cpus;

fn determine_buffer_size() -> usize {
    const DEFAULT_SIZE: usize = 5 * 1024 * 1024; // Default to 1MB
    const MAX_SIZE: usize = 10 * 1024 * 1024; // Cap at 10MB
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
    #[clap(short, long)]
    csv: Option<String>,

    /// Directory to be processed. Non-recursively.
    #[clap(short, long)]
    directory: Option<String>,

    /// Numbers of threads to be used, otherwise the program will decide.
    /// 
    /// It will decide based on the number of available logical threads, physical cpus, and number of files 
    /// to be processed.
    #[clap(short, long)]
    threads: Option<usize>,

    /// Legacy output
    #[clap(short, long)]
    legacy: bool,

    /// FASTA FILE[s] to be processed [wildcards would work here].
    /// 
    /// Inside a zip file, only .fa .fasta .fna files will be processed.
    /// gz files are assumed to be compressed using gzip, wrong results can come out of a gz file compressed using bzip
    #[clap(name = "FASTA FILE")]
    files: Vec<String>,
}

#[derive(Default, Clone, Debug)] // Added Debug derive
struct AnalysisResults {
    filename: String,
    total_length: usize,
    sequence_count: usize,
    gc_count: usize,
    n_count: usize,
    n25: usize,
    n25_sequence_count: usize,
    n50: usize,
    n50_sequence_count: usize,
    n75: usize,
    n75_sequence_count: usize,
    largest_contig: usize,
    shortest_contig: usize,
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

const VALID_FILES: [&str; 3] = ["fa", "fasta", "fna"];

fn get_fasta_files_from_directory(dir: &str) -> std::io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file() {
            if let Some(ext) = path.extension().and_then(|ext| ext.to_str()) {
                if VALID_FILES.contains(&ext) {
                    files.push(path);
                } else if ["gz", "xz", "bz2", "bgz", "bgzip", "zip"].contains(&ext) {
                    files.push(path);
                }
            }
        }
    }
    Ok(files)
}

fn process_files(files: Vec<PathBuf>, threads: Option<usize>) -> Vec<AnalysisResults> {
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
                    "gz" => process_gz_file(file, buffer_size).unwrap_or_default(),
                    "zip" => process_zip_file(file, buffer_size).unwrap_or_default(),
                    "xz" => process_xz_file(file, buffer_size).unwrap_or_default(),
                    "bz2" => process_bz2_file(file, buffer_size).unwrap_or_default(),
                    "bgz" | "bgzip" => process_bgzip_file(file, buffer_size).unwrap_or_default(),
                    _ if VALID_FILES.contains(&ext) => {
                        process_fasta_file(file, buffer_size).unwrap_or_default()
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

fn process_fasta_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let file = File::open(file)?;
    let reader = BufReader::with_capacity(buffer_size, file);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

fn process_gz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let file = File::open(file)?;
    let gz = GzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, gz);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

fn process_zip_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let file = File::open(file)?;
    let buf_reader = BufReader::with_capacity(buffer_size, file);
    let mut archive = ZipArchive::new(buf_reader)?;
    let mut all_results = Vec::new();

    for i in 0..archive.len() {
        let zip_file = archive.by_index(i)?;
        if zip_file.is_file() {
            let file_name = zip_file.name().to_owned();
            if VALID_FILES.iter().any(|&ext| file_name.ends_with(ext)) {
                let mut result = AnalysisResults {
                    filename: Path::new(&file_name)
                        .file_name()
                        .unwrap()
                        .to_string_lossy()
                        .to_string(),
                    shortest_contig: usize::MAX,
                    ..Default::default()
                };
                let reader = BufReader::with_capacity(buffer_size, zip_file);
                if let Err(e) = process_reader(reader, &mut result) {
                    eprintln!("Error processing {}: {}", file_name, e);
                    continue; // Skip this file but continue processing others
                };
                all_results.push(result);
            }
        }
    }

    Ok(all_results)
}

fn process_reader<R: Read>(
    mut reader: BufReader<R>,
    results: &mut AnalysisResults,
) -> std::io::Result<()> {
    let mut lengths = Vec::with_capacity(250);
    let mut current_sequence_length = 0;
    let mut line = Vec::with_capacity(128);
    let offset;

    if reader.read_until(b'\n', &mut line)? > 0 {
        results.sequence_count += 1;
        //Assuming the first line is a header line and starts with '>'
        offset = {
            if line.ends_with(b"\r\n") {
                Some(2) // Exclude the newline character
            } else if line.ends_with(b"\n") {
                Some(1) // Exclude the newline characters
            } else {
                None // No newline characters?
            }
        };
        line.clear();
    } else {
        return Ok(()); // Nothing read from the file
    };
    let Some(offset) = offset else {
        return Ok(()); // No newline characters?
    };

    while reader.read_until(b'\n', &mut line)? > 0 {
        // Already processed the first line
        if line.first() == Some(&b'>') {
            results.total_length += current_sequence_length;
            results.largest_contig = results.largest_contig.max(current_sequence_length);
            results.shortest_contig = results.shortest_contig.min(current_sequence_length);
            lengths.push(current_sequence_length);
            current_sequence_length = 0;
            results.sequence_count += 1;
        } else {
            current_sequence_length += process_sequence_line(&line, results, offset);
        }
        line.clear();
    }

    if current_sequence_length > 0 {
        results.total_length += current_sequence_length;
        results.largest_contig = results.largest_contig.max(current_sequence_length);
        results.shortest_contig = results.shortest_contig.min(current_sequence_length);
        lengths.push(current_sequence_length);
    }

    calc_nq_stats(&lengths, results);
    Ok(())
}

fn process_sequence_line(line: &[u8], results: &mut AnalysisResults, offset: usize) -> usize {
    results.gc_count += bytecount::count(line, b'G')
        + bytecount::count(line, b'g')
        + bytecount::count(line, b'C')
        + bytecount::count(line, b'c');
    results.n_count += bytecount::count(line, b'N') + bytecount::count(line, b'n');
    if line.ends_with(b"\n") {
        line.len() - offset // Exclude the newline character
    } else {
        line.len()
    }
}

fn calc_nq_stats(lengths: &[usize], results: &mut AnalysisResults) {
    let total_length: usize = lengths.iter().sum();
    let mut cumulative_length = 0;
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    for (i, &length) in sorted_lengths.iter().enumerate() {
        cumulative_length += length;
        if results.n25 == 0 && cumulative_length >= total_length / 4 {
            results.n25 = length;
            results.n25_sequence_count = i + 1;
        }
        if results.n50 == 0 && cumulative_length >= total_length / 2 {
            results.n50 = length;
            results.n50_sequence_count = i + 1;
        }
        if results.n75 == 0 && cumulative_length >= total_length * 3 / 4 {
            results.n75 = length;
            results.n75_sequence_count = i + 1;
            break;
        }
    }
}

fn print_results(results: &AnalysisResults, legacy: bool) {
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

fn append_to_csv(results: &[AnalysisResults], csv_filename: &str) -> io::Result<()> {
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

fn process_xz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let file = File::open(file)?;
    let xz = XzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, xz);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

fn process_bz2_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let file = File::open(file)?;
    let bz = BzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, bz);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

fn process_bgzip_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let file = File::open(file)?;
    let mut reader = bgzf::Reader::new(file);
    let mut buffer = Vec::new();

    // Read entire decompressed content
    reader.read_to_end(&mut buffer)?;

    let reader = BufReader::with_capacity(buffer_size, &buffer[..]);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn it_works() {
        let mut files_to_process = Vec::new();

        if let Ok(files) = get_fasta_files_from_directory(&"./test/") {
            files_to_process.extend(files);
        }

        let results = process_files(files_to_process, None);

        let csv_file = "test/attempt.csv";

        append_to_csv(&results, &csv_file).expect("Failed to write CSV");
        let thing = fs::read_to_string("test/test.csv").unwrap();
        let compare = fs::read_to_string(csv_file).unwrap();
        assert_eq!(thing, compare);
        let _ = fs::remove_file(csv_file);
    }
}
