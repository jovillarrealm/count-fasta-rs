use clap::Parser;
use flate2::read::GzDecoder;
use futures::stream::{self, StreamExt};
use memchr::{memchr2_iter, memchr_iter, memrchr2};
use parking_lot::Mutex;
use std::env;
use std::io::{BufRead, BufReader as StdBufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{self, AsyncBufReadExt, AsyncWriteExt, BufReader as TokioBufReader};
use zip::read::ZipArchive;

fn determine_buffer_size() -> usize {
    let default_size = 8 * 1024; // Default to 1kB
    let max_size = 10 * 1024 * 1024; // Cap at 10MB
    match env::var("BUFFER_SIZE") {
        Ok(val) => val.parse().unwrap_or(default_size).min(max_size),
        Err(_) => default_size,
    }
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

#[tokio::main]
async fn main() -> io::Result<()> {
    let args = Args::parse();

    let mut files_to_process = Vec::new();

    if let Some(dir) = args.directory {
        files_to_process.extend(get_fasta_files_from_directory(&dir).await?);
    }
    files_to_process.extend(args.files.into_iter().map(PathBuf::from));

    let results = process_files(files_to_process).await?;

    if let Some(csv_file) = args.csv {
        append_to_csv(&results, &csv_file).await?;
    } else {
        for result in results {
            print_results(&result);
        }
    }

    Ok(())
}

async fn get_fasta_files_from_directory(dir: &str) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    let mut entries = tokio::fs::read_dir(dir).await?;
    while let Some(entry) = entries.next_entry().await? {
        let path = entry.path();
        if path.is_file() {
            let filename = path.file_name().and_then(|file_name| file_name.to_str());
            if let Some(filename) = filename {
                if filename.ends_with(".fna.gz")
                    || filename.ends_with(".fa.gz")
                    || filename.ends_with(".fasta.gz")
                    || filename.ends_with(".zip")
                    || filename.ends_with(".fna")
                    || filename.ends_with(".fa")
                    || filename.ends_with(".fasta")
                {
                    files.push(path);
                }
            }
        }
    }
    Ok(files)
}

async fn process_files(files: Vec<PathBuf>) -> io::Result<Vec<AnalysisResults>> {
    let results = Arc::new(Mutex::new(Vec::with_capacity(50_000)));
    let buffer_size = determine_buffer_size();

    stream::iter(files)
        .map(|file| {
            let results = Arc::clone(&results);
            async move {
                let file_str = file.to_str().unwrap();
                let mut local_result = AnalysisResults {
                    filename: Path::new(&file_str)
                        .file_name()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .to_string(),
                    shortest_contig: usize::MAX,
                    ..Default::default()
                };

                let result: io::Result<()> = match file_str {
                    f if f.ends_with(".gz") => process_gz_file(&file, &mut local_result).await,
                    f if f.ends_with(".zip") => process_zip_file(&file, &mut local_result).await,
                    _ => process_fasta_file(&file, &mut local_result, buffer_size).await,
                };

                if let Err(e) = result {
                    eprintln!("Error processing file {}: {}", file_str, e);
                } else {
                    results.lock().push(local_result);
                }
            }
        })
        .buffer_unordered(
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
        )
        .collect::<Vec<_>>()
        .await;

    Ok(Arc::try_unwrap(results).unwrap().into_inner())
}

async fn process_fasta_file(
    file: &Path,
    results: &mut AnalysisResults,
    buffer_size: usize,
) -> io::Result<()> {
    let file = File::open(file).await?;
    let mut reader = TokioBufReader::with_capacity(buffer_size, file);
    let mut buffer = Vec::with_capacity(buffer_size);
    let mut lengths = Vec::with_capacity(500);
    let mut current_sequence_length = 0;
    let mut delimiter = b'\n';
    let mut sequence_length ;

    while reader.read_until(delimiter, &mut buffer).await? > 0 {
        if delimiter == b'\n' {
            if current_sequence_length > 0 {
                results.total_length += current_sequence_length;
                results.largest_contig = results.largest_contig.max(current_sequence_length);
                results.shortest_contig = results.shortest_contig.min(current_sequence_length);
                lengths.push(current_sequence_length);
                current_sequence_length = 0;
            }
            results.sequence_count += 1;
            buffer.clear();
            delimiter = b'>';
        } else {
            sequence_length = process_sequence_bytes(&buffer, results);
            current_sequence_length += sequence_length;
            delimiter = b'\n';
        }
    }

    if current_sequence_length > 0 {
        results.total_length += current_sequence_length;
        results.largest_contig = results.largest_contig.max(current_sequence_length);
        results.shortest_contig = results.shortest_contig.min(current_sequence_length);
        lengths.push(current_sequence_length);
        current_sequence_length = 0;
    }

    if current_sequence_length > 0 {
        lengths.push(current_sequence_length);
    }

    calc_nq_stats(&mut lengths, results);
    Ok(())
}

async fn process_gz_file(file: &Path, results: &mut AnalysisResults) -> io::Result<()> {
    let file = std::fs::File::open(file)?;
    let gz = GzDecoder::new(file);
    let reader = StdBufReader::new(gz);
    process_reader(reader, results)
}

async fn process_zip_file(file: &Path, results: &mut AnalysisResults) -> io::Result<()> {
    let file = std::fs::File::open(file)?;
    let mut archive = ZipArchive::new(file)?;
    // Iterate through each file in the ZIP archive
    for i in 0..archive.len() {
        let file = archive.by_index(i)?;
        if file.is_file() {
            let file_name = file.name();

            if file_name.ends_with(".fasta")
                || file_name.ends_with(".fa")
                || file_name.ends_with(".fna")
            {
                // Process the file (example: pass to process_fasta_file)
                let reader = StdBufReader::new(file);
                return process_reader(reader, results);
            }
        }
    }
    Ok(())
}

fn process_reader<R: Read>(
    reader: StdBufReader<R>,
    results: &mut AnalysisResults,
) -> io::Result<()> {
    let mut lengths = Vec::with_capacity(250);
    let mut current_sequence_length = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('>') {
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
    }

    if current_sequence_length > 0 {
        results.total_length += current_sequence_length;
        results.largest_contig = results.largest_contig.max(current_sequence_length);
        results.shortest_contig = results.shortest_contig.min(current_sequence_length);
        lengths.push(current_sequence_length);
    }

    calc_nq_stats(&mut lengths, results);
    Ok(())
}

fn process_sequence_bytes(buffer: &[u8], results: &mut AnalysisResults) -> usize {
    results.gc_count += memchr2_iter(b'G', b'C', buffer).count() as usize;
    results.n_count += memchr_iter(b'N', buffer).count() as usize;
    let noise = memchr2_iter(b'>', b'\n', buffer).count();
    buffer.len() - noise
}

fn process_sequence_line(line: &str, results: &mut AnalysisResults) -> usize {
    let line = line.trim();
    let line_bytes = line.as_bytes();
    results.gc_count += memchr2_iter(b'G', b'C', line_bytes).count() as usize;
    results.n_count += memchr_iter(b'N', line_bytes).count() as usize;
    line_bytes.len()
}

fn calc_nq_stats(lengths: &mut [usize], results: &mut AnalysisResults) {
    let total_length: usize = lengths.iter().sum();
    let mut cumulative_length = 0;
    lengths.sort_unstable_by(|a, b| b.cmp(a));

    for (i, &length) in lengths.iter().enumerate() {
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

async fn append_to_csv(results: &[AnalysisResults], csv_filename: &str) -> io::Result<()> {
    let csv_exists = Path::new(csv_filename).exists();
    let mut file = tokio::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(csv_filename)
        .await?;

    if !csv_exists {
        let header = "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n";
        file.write_all(header.as_bytes()).await?;
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
            file.write_all(buffer.as_bytes()).await?;
            buffer.clear();
        }
    }

    if !buffer.is_empty() {
        file.write_all(buffer.as_bytes()).await?;
    }

    file.flush().await?;
    Ok(())
}
