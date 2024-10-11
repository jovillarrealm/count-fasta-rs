//! # FASTA Sequence Analyzer
//!
//! This program analyzes one or more FASTA format files containing DNA sequences.
//! It calculates various statistics such as sequence length, GC content,
//! and N50 values. Results can be displayed on the console and optionally
//! appended to a CSV file for further analysis.

//use std::collections::HashMap;
use seq_io::fasta::{Reader, RefRecord};
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{self, stdout, BufReader, BufWriter, Write};
use std::os::unix::io::AsRawFd;
use std::path::Path;
use flate2::read::GzDecoder;
use zip::read::ZipArchive;

/// Holds the results of the sequence analysis.
#[derive(Default)]
struct AnalysisResults<'a> {
    /// Names of the analyzed FASTA files.
    filename: &'a str,
    /// Total length of all sequences combined.
    total_length: usize,
    /// Number of sequences in all files.
    sequence_count: usize,
    /// Total count of G and C bases.
    gc_count: usize,
    /// Total count of N bases.
    n_count: usize,
    /// N25 statistic (length at 25% of total sequence length).
    n25: usize,
    n25_sequence_count: usize,
    /// N50 statistic (length at 50% of total sequence length).
    n50: usize,
    n50_sequence_count: usize,
    /// N75 statistic (length at 75% of total sequence length).
    n75: usize,
    n75_sequence_count: usize,
    /// Length of the largest contig (sequence).
    largest_contig: usize,
    /// Length of the shortest contig (sequence).
    shortest_contig: usize,
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} [-c csv_file] <fasta_file1>", args[0]);
        std::process::exit(1);
    }

    let (csv_filename, fasta_file) = parse_args(&args);

    let file = File::open(fasta_file)?;
    let buffer_size = 1024*1024*16;
    if fasta_file.ends_with(".gz") {
        let gz_reader = GzDecoder::new(file);
        let buf_reader = BufReader::with_capacity(buffer_size, gz_reader);
        process_file(buf_reader, fasta_file, csv_filename)?;

    } else if fasta_file.ends_with(".zip") {
        let mut archive = ZipArchive::new(file)?;
        let file_in_zip = archive.by_index(0)?; // Assuming the first file
        let buf_reader = BufReader::with_capacity(buffer_size, file_in_zip);
        process_file(buf_reader, fasta_file, csv_filename)?;
    } else {
        let buf_reader = BufReader::with_capacity(buffer_size, file);
        process_file(buf_reader, fasta_file, csv_filename)?;
    }
    Ok(())
}

fn process_file<R: io::Read>(buf_reader: BufReader<R>, fasta_file: &str, csv_filename: Option<&str>) -> Result<(), io::Error> {
    let mut fasta_reader = Reader::new(buf_reader);
    let mut results = AnalysisResults {
        filename: fasta_file,
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let mut lengths = Vec::with_capacity(500);
    while let Some(result) = fasta_reader.next() {
        let record = result.expect("Failed to get RefRecord of secuence") ;
        process_sequence(&record, &mut results, &mut lengths)?;
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    if let Some(csv_file) = csv_filename {
        calc_n50_stats(& lengths, &mut results);
        append_to_csv(&results, csv_file)?;
    } else {
        calc_nq_stats(&lengths, &mut results);
        print_results(&results);
    }
    Ok(())
}
/// Analyzes the given sequences and computes various statistics.
///
/// # Arguments
///
/// * `sequences` - A slice of `Sequence` structs to analyze
/// * `interval_size` - The interval size for the length histogram
///
/// # Returns
///
/// An `AnalysisResults` struct containing the computed statistics.
///

fn process_sequence(record: &RefRecord, results: &mut AnalysisResults, lengths: &mut Vec<usize>) -> io::Result<()> {
    results.sequence_count += 1;
    let mut sequence_length:usize = 0;
    let mut gc_count = 0;
    let mut n_count = 0;
    
    for line in record.seq_lines() {
        sequence_length += line.len();

        for &c in line {
            match c {
                b'G' | b'C' | b'g' | b'c' => gc_count += 1,
                b'N' | b'n' => n_count += 1,
                _ => (),
            }
        }
    }

    results.total_length += sequence_length;
    results.gc_count += gc_count;
    results.n_count += n_count;
    results.largest_contig = results.largest_contig.max(sequence_length);
    results.shortest_contig = results.shortest_contig.min(sequence_length);
    
    lengths.push(sequence_length);
    Ok(())
}

/// Calculate N25, N50, N75 stats
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `interval_size` - The interval size used for the length histogram
fn calc_nq_stats(lengths: &[usize], results: &mut AnalysisResults) {

    let total_length = results.total_length;
    let mut cumulative_length = 0;
    let mut cumulative_sequences = 0;

    for &length in lengths {
        cumulative_length += length;
        cumulative_sequences += 1;

        if results.n25 == 0 && cumulative_length >= total_length / 4 {
            results.n25 = length;
            results.n25_sequence_count = cumulative_sequences;
        }
        if results.n50 == 0 && cumulative_length >= total_length / 2 {
            results.n50 = length;
            results.n50_sequence_count = cumulative_sequences;
        }
        if results.n75 == 0 && cumulative_length >= total_length * 3 / 4 {
            results.n75 = length;
            results.n75_sequence_count = cumulative_sequences;
            break;
        }
    }
}


fn calc_n50_stats(lengths: &[usize], results: &mut AnalysisResults) {
    let total_length = results.total_length;
    let mut cumulative_length = 0;
    let mut cumulative_sequences = 0;

    for &length in lengths {
        cumulative_length += length;
        cumulative_sequences += 1;

        if results.n25 == 0 && cumulative_length >= total_length / 4 {
            results.n25 = length;
            results.n25_sequence_count = cumulative_sequences;
        }
        if results.n50 == 0 && cumulative_length >= total_length / 2 {
            results.n50 = length;
            results.n50_sequence_count = cumulative_sequences;
            break;
        }
        if results.n75 == 0 && cumulative_length >= total_length * 3 / 4 {
            results.n75 = length;
            results.n75_sequence_count = cumulative_sequences;
        }
    }
}

/// Prints the analysis results to the console.
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `interval_size` - The interval size used for the length histogram
//fn print_results(results: &AnalysisResults) {
fn print_results(results: &AnalysisResults) {
    let mut stdout = BufWriter::new(stdout());
    writeln!(
        &mut stdout,
        "\nTotal length of sequence:\t{} bp",
        results.total_length
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Total number of sequences:\t{}",
        results.sequence_count
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Average contig length is:\t{} bp",
        results.total_length / results.sequence_count
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Largest contig:\t\t{} bp",
        results.largest_contig
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Shortest contig:\t\t{} bp",
        results.shortest_contig
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N25 stats:\t\t\t25% of total sequence length is contained in the {} sequences >= {} bp",
        results.n25_sequence_count, results.n25
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N50 stats:\t\t\t50% of total sequence length is contained in the {} sequences >= {} bp",
        results.n50_sequence_count, results.n50
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N75 stats:\t\t\t75% of total sequence length is contained in the {} sequences >= {} bp",
        results.n75_sequence_count, results.n75
    )
    .unwrap();
    writeln!(&mut stdout, "Total GC count:\t\t\t{} bp", results.gc_count).unwrap();
    writeln!(
        &mut stdout,
        "GC %:\t\t\t\t{:.2} %",
        (results.gc_count as f64 / results.total_length as f64) * 100.0
    )
    .unwrap();
    writeln!(&mut stdout, "Number of Ns:\t\t\t{}", results.n_count).unwrap();
    writeln!(
        &mut stdout,
        "Ns %:\t\t\t\t{:.2} %",
        (results.n_count as f64 / results.total_length as f64) * 100.0
    )
    .unwrap();
    stdout.flush().unwrap();
}
/// Appends the analysis results to a CSV file.
///
/// If the file doesn't exist, it will be created and a header row will be written.
/// If the file exists, the results will be appended as a new row.
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `csv_filename` - The name of the CSV file to append to
///
/// # Returns
///
/// An `io::Result<()>` indicating success or containing an error if the operation failed.
fn append_to_csv(results: &AnalysisResults, csv_filename: &str) -> io::Result<()> {
    let file_exists = Path::new(csv_filename).exists();
    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(csv_filename)?;
    
    // Use `flock` to lock the file before writing
    let fd = file.as_raw_fd();
    unsafe { libc::flock(fd, libc::LOCK_EX) }; // Acquire exclusive lock


    let mut csv_writer = CsvWriter::new(BufWriter::new(file));

    if !file_exists {
        csv_writer.write_header()?;
    }
    csv_writer.write_record(results)?;

    // Release the file lock
    unsafe { libc::flock(fd, libc::LOCK_UN) }; // Unlock the file
    Ok(())
}

fn parse_args(args: &[String]) -> (Option<&str>, &str) {
    let mut csv_filename = None;
    let mut fasta_file = "";
    let mut i = 1;

    while i < args.len() {
        match args[i].as_str() {
            "-c" | "--csv" if i + 1 < args.len() => {
                csv_filename = Some(&args[i + 1]);
                i += 2;
            }
            _ => {
                if fasta_file.is_empty() {
                    fasta_file = &args[i];
                }
                i += 1;
            }
        }
    }

    (csv_filename.map(|x| x.as_str()), fasta_file)
}



struct CsvWriter<W: Write> {
    writer: W,
}

impl<W: Write> CsvWriter<W> {
    fn new(writer: W) -> Self {
        CsvWriter { writer }
    }

    fn write_header(&mut self) -> io::Result<()> {
        self.writer.write_all(b"filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n")
    }

    fn write_record(&mut self, results: &AnalysisResults) -> io::Result<()> {
        let filename = Path::new(&results.filename).file_name().unwrap().to_str().unwrap();
        write!(self.writer, "{};", filename)?;
        write!(self.writer, "{};", results.total_length)?;
        write!(self.writer, "{};", results.sequence_count)?;
        write!(self.writer, "{};", results.total_length as f64 / results.sequence_count as f64)?;
        write!(self.writer, "{};", results.largest_contig)?;
        write!(self.writer, "{};", results.shortest_contig)?;
        write!(self.writer, "{};", results.n50)?;
        write!(self.writer, "{};", (results.gc_count as f64 / results.total_length as f64) * 100.0)?;
        write!(self.writer, "{};", results.n_count)?;
        writeln!(self.writer, "{}", (results.n_count as f64 / results.total_length as f64) * 100.0)?;
        Ok(())
    }
}
