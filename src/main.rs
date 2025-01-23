use clap::Parser;
use crossbeam_channel::unbounded;
use flate2::read::GzDecoder;
use memchr;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::BufRead;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use zip::ZipArchive;

// Increased default buffer size for HDD optimization
fn determine_buffer_size() -> usize {
    env::var("BUFFER_SIZE")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1024 * 1024) // Default 1MB
        .min(10 * 1024 * 1024) // Max 10MB
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long)]
    csv: Option<String>,

    #[clap(short, long)]
    directory: Option<String>,

    #[clap(name = "FASTA FILE (.fa .fasta .fna .zip .gz)")]
    files: Vec<String>,
}

#[derive(Default, Clone, Debug)]
struct AnalysisResults {
    filename: String,
    total_length: usize,
    sequence_count: usize,
    gc_count: usize,
    n_count: usize,
    lengths: Vec<usize>,
    n25: usize,
    n25_sequence_count: usize,
    n50: usize,
    n50_sequence_count: usize,
    n75: usize,
    n75_sequence_count: usize,
    largest_contig: usize,
    shortest_contig: usize,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let files = collect_files(&args)?;
    let (sender, receiver) = unbounded();

    // Configure thread pool for HDD-friendly parallelism
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus::get().min(6)) // Limit threads for HDD
        .build()
        .unwrap()
        .install(|| {
            files.into_par_iter().for_each_with(sender, |s, path| {
                let mut result = match process_file(&path) {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Error processing {}: {}", path.display(), e);
                        AnalysisResults::default()
                    }
                };
                result.filename = path.file_name().unwrap().to_string_lossy().into_owned();
                let _ = s.send(result);
            });
        });

    let results: Vec<_> = receiver.try_iter().collect();

    if let Some(csv_path) = args.csv {
        write_csv(&results, &csv_path)?;
    } else {
        results.iter().for_each(print_results);
    }

    Ok(())
}

fn collect_files(args: &Args) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    if let Some(dir) = &args.directory {
        for entry in std::fs::read_dir(dir)? {
            let path = entry?.path();
            if path.is_file() && is_fasta_file(&path) {
                files.push(path);
            }
        }
    }

    files.extend(args.files.iter().map(PathBuf::from));
    Ok(files)
}

fn is_fasta_file(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ["fa", "fasta", "fna", "gz", "zip"].contains(&ext))
        .unwrap_or(false)
}

fn process_file(path: &Path) -> io::Result<AnalysisResults> {
    let mut results = AnalysisResults {
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let buffer_size = determine_buffer_size();

    match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => {
            let file = File::open(path)?;
            process_fasta_chunked(GzDecoder::new(file), &mut results, buffer_size)?;
        }
        Some("zip") => {
            let file = File::open(path)?;
            let mut archive = ZipArchive::new(file)?;
            // Process first FASTA file found in archive
            for i in 0..archive.len() {
                let file = archive.by_index(i)?;
                if is_fasta_file(Path::new(file.name())) {
                    process_fasta_chunked(file, &mut results, buffer_size)?;
                    break;
                }
            }
        }
        _ => {
            let file = File::open(path)?;
            process_fasta_chunked(file, &mut results, buffer_size)?;
        }
    }

    calc_nq_stats(&mut results);
    Ok(results)
}

fn process_fasta_chunked<R: Read>(
    reader: R,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> io::Result<()> {
    let mut reader = BufReader::with_capacity(buffer_size, reader);
    let mut buffer = Vec::with_capacity(buffer_size);
    let mut carry_over = Vec::new();
    let mut in_sequence = false;

    loop {
        buffer.clear();
        let bytes_read = reader.read_until(b'>', &mut buffer)?;
        if bytes_read == 0 {
            break;
        }

        // Combine carry_over with new buffer contents
        let mut combined = carry_over.clone();
        combined.extend_from_slice(&buffer);

        let headers = memchr::memchr_iter(b'>', &combined)
            .filter(|&pos| pos == 0 || combined[pos - 1] == b'\n');

        let mut last_pos = 0;
        let mut processed_headers = 0;
        let headers2 = headers.clone();

        for header_pos in headers2 {
            if in_sequence {
                // Process sequence data between last header and this one
                process_sequence(&combined[last_pos..header_pos], results);
            }

            // Skip the entire header line
            let header_end = combined[header_pos..]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| header_pos + p + 1)
                .unwrap_or(combined.len());

            last_pos = header_end;
            in_sequence = true;
            processed_headers += 1;
        }

        // Process remaining data after last header
        if in_sequence && last_pos < combined.len() {
            process_sequence(&combined[last_pos..], results);
        }

        // Update carry_over with unprocessed data
        carry_over.clear();
        if !headers.last().is_some_and(|hp| hp >= combined.len()) {
            carry_over.extend_from_slice(&combined);
        }
    }

    Ok(())
}

fn process_sequence(data: &[u8], results: &mut AnalysisResults) {
    let mut length = 0;
    let mut gc = 0;
    let mut n = 0;

    for &byte in data {
        match byte {
            b'>' => break,             // Stop processing if we encounter a new header
            b'\n' | b'\r' => continue, // Skip newline characters
            _ => {
                length += 1;
                match byte.to_ascii_uppercase() {
                    b'G' | b'C' => gc += 1,
                    b'N' => n += 1,
                    _ => (),
                }
            }
        }
    }

    if length > 0 {
        results.sequence_count += 1;
        results.total_length += length;
        results.gc_count += gc;
        results.n_count += n;
        results.lengths.push(length);
        results.largest_contig = results.largest_contig.max(length);
        results.shortest_contig = results.shortest_contig.min(length);
    }
}
fn calc_nq_stats(results: &mut AnalysisResults) {
    let total_length = results.total_length;
    let mut sorted_lengths = results.lengths.clone();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    let mut cumulative = 0;
    for (i, &len) in sorted_lengths.iter().enumerate() {
        cumulative += len;

        if results.n25 == 0 && cumulative >= total_length / 4 {
            results.n25 = len;
            results.n25_sequence_count = i + 1;
        }
        if results.n50 == 0 && cumulative >= total_length / 2 {
            results.n50 = len;
            results.n50_sequence_count = i + 1;
        }
        if cumulative >= total_length * 3 / 4 {
            results.n75 = len;
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

fn write_csv(results: &[AnalysisResults], path: &str) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "filename;total_length;sequences;largest_contig;shortest_contig;n50;gc_percent;n_percent"
    )?;

    for result in results {
        let gc_pct = (result.gc_count as f64 / result.total_length as f64) * 100.0;
        let n_pct = (result.n_count as f64 / result.total_length as f64) * 100.0;
        writeln!(
            writer,
            "{};{};{};{};{};{};{:.6};{:.6}",
            result.filename,
            result.total_length,
            result.sequence_count,
            result.largest_contig,
            result.shortest_contig,
            result.n50,
            gc_pct,
            n_pct
        )?;
    }

    writer.flush()
}
