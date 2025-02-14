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
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long)]
    csv: Option<String>,

    #[clap(short, long)]
    directory: Option<String>,

    #[clap(
        name = "FASTA FILE (.fa .fasta .fna .zip(multi file archive) .gz .xz .bz2 .bgz .bgzip)"
    )]
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

    let results = process_files(files_to_process);

    if let Some(csv_file) = args.csv {
        append_to_csv(&results, &csv_file).expect("Failed to write CSV");
    } else {
        for result in results {
            print_results(&result);
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
                    path.clone().file_name().and_then(|osf_name| {
                        osf_name.to_str().and_then(|f_name| {
                            VALID_FILES
                                .iter()
                                .any(|&v| f_name.contains(v))
                                .then(|| files.push(path))
                        })
                    });
                }
            }
        }
    }
    Ok(files)
}

fn process_files(files: Vec<PathBuf>) -> Vec<AnalysisResults> {
    let buffer_size = determine_buffer_size();
    let usable_threads = (num_cpus::get() as f32 * 0.9).round() as usize;
    let available_threads = max(
        1,
        min(usable_threads, files.len()),
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(available_threads)
        .stack_size(1024 * 1024)
        .build()
        .unwrap();
    pool.install(|| {
        files
            .par_iter()
            .filter_map(|file| {
                let mut local_result = AnalysisResults {
                    filename: file.file_name().unwrap().to_string_lossy().to_string(),
                    shortest_contig: usize::MAX,
                    ..Default::default()
                };

                let result = match file.extension()?.to_str()? {
                    "gz" => process_gz_file(file, &mut local_result, buffer_size),
                    "zip" => process_zip_file(file, &mut local_result, buffer_size),
                    "xz" => process_xz_file(file, &mut local_result, buffer_size),
                    "bz2" => process_bz2_file(file, &mut local_result, buffer_size),
                    "bgz" | "bgzip" => process_bgzip_file(file, &mut local_result, buffer_size),
                    ext if VALID_FILES.contains(&ext) => {
                        process_fasta_file(file, &mut local_result, buffer_size)
                    }
                    _ => Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "Unsupported file type",
                    )),
                };

                result.ok().map(|_| local_result)
            })
            .collect()
    })
}

fn process_fasta_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let reader = BufReader::with_capacity(buffer_size, file);
    process_reader(reader, results)
}

fn process_gz_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let gz = GzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, gz);
    process_reader(reader, results)
}

fn process_zip_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let buf_reader = BufReader::with_capacity(buffer_size, file);
    let mut archive = ZipArchive::new(buf_reader)?;
    for i in 0..archive.len() {
        let mut zip_file = archive.by_index(i)?;
        if zip_file.is_file() {
            let file_name = zip_file.name();
            if VALID_FILES.iter().any(|&ext| file_name.ends_with(ext)) {
                let reader = BufReader::with_capacity(buffer_size, &mut zip_file);
                return process_reader(reader, results);
            }
        }
    }
    Ok(())
}

fn process_reader<R: Read>(
    mut reader: BufReader<R>,
    results: &mut AnalysisResults,
) -> std::io::Result<()> {
    let mut lengths = Vec::with_capacity(250);
    let mut current_sequence_length = 0;
    let mut line = Vec::with_capacity(128);

    while reader.read_until(b'\n', &mut line)? > 0 {
        if line.first() == Some(&b'>') {
            if current_sequence_length > 0 {
                results.total_length += current_sequence_length;
                results.largest_contig = results.largest_contig.max(current_sequence_length);
                results.shortest_contig = results.shortest_contig.min(current_sequence_length);
                lengths.push(current_sequence_length);
                current_sequence_length = 0;
            }
            results.sequence_count += 1;
        } else {
            current_sequence_length += process_sequence_line(&line, results);
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

fn process_sequence_line(line: &[u8], results: &mut AnalysisResults) -> usize {
    results.gc_count += bytecount::count(line, b'G')
        + bytecount::count(line, b'g')
        + bytecount::count(line, b'C')
        + bytecount::count(line, b'c');
    results.n_count += bytecount::count(line, b'N') + bytecount::count(line, b'n');
    if line.ends_with(b"\n") {
        line.len() - 1 // Exclude the newline character
    } else if line.ends_with(b"\r\n") {
        line.len() - 2 // Exclude the newline characters
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

fn print_results(results: &AnalysisResults) {
    println!("\nTotal length of sequence:\t{} bp", results.total_length);
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
            "{};{};{};{:.7};{};{};{};{:.7};{};{:.7}\n",
            result.filename,
            result.total_length,
            result.sequence_count,
            result.total_length as f64 / result.sequence_count as f64,
            result.largest_contig,
            result.shortest_contig,
            result.n50,
            (result.gc_count as f64 / result.total_length as f64) * 100.0,
            result.n_count,
            (result.n_count as f64 / result.total_length as f64) * 100.0,
        );
        buffer.push_str(&line);

        // Write in chunks to avoid holding too much in memory
        if buffer.len() > 64 * 1024 {
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

fn process_xz_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let xz = XzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, xz);
    process_reader(reader, results)
}

fn process_bz2_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let bz = BzDecoder::new(file);
    let reader = BufReader::with_capacity(buffer_size, bz);
    process_reader(reader, results)
}

fn process_bgzip_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> std::io::Result<()> {
    let file = File::open(file)?;
    let mut reader = bgzf::Reader::new(file);
    let mut buffer = Vec::new();

    // Read entire decompressed content
    reader.read_to_end(&mut buffer)?;

    let reader = BufReader::with_capacity(buffer_size, &buffer[..]);
    process_reader(reader, results)
}
