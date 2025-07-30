// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use liblzma::read::XzDecoder;
use noodles_bgzf as bgzf;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use zip::read::ZipArchive;

pub const VALID_FILES: [&str; 3] = ["fa", "fasta", "fna"];
pub const VALID_COMPRESSION: [&str; 7] = ["gz", "xz", "bz2", "bgz", "bgzip", "zip", "naf"];

#[derive(Default, Clone, Debug)]
pub struct AnalysisResults {
    pub filename: String,
    pub total_length: usize,
    pub sequence_count: usize,
    pub gc_count: usize,
    pub n_count: usize,
    pub n25: usize,
    pub n25_sequence_count: usize,
    pub n50: usize,
    pub n50_sequence_count: usize,
    pub n75: usize,
    pub n75_sequence_count: usize,
    pub largest_contig: usize,
    pub shortest_contig: usize,
}

pub fn process_xz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
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

pub fn process_bz2_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
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

pub fn process_bgzip_file(
    file: &Path,
    buffer_size: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let mut reader = File::open(file).map(bgzf::io::Reader::new)?;
    let mut buffer = Vec::new();

    // Read entire decompressed content
    reader.read_to_end(&mut buffer)?;

    let reader = BufReader::with_capacity(buffer_size, &buffer[..]);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

pub fn process_fasta_file(
    file: &Path,
    buffer_size: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
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

pub fn process_naf_file(file: &Path) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults {
        filename: file.file_name().unwrap().to_string_lossy().to_string(),
        shortest_contig: usize::MAX,
        ..Default::default()
    };
    let decoder = nafcodec::Decoder::from_path(file)
        .expect("failed to open nucleotide archive");

    // Process naf file
    let mut lengths = Vec::with_capacity(250);
    let offset =0;

    for may_seq in decoder {
        let seq = may_seq.unwrap_or_else(|_| panic!("{file:?} had bad data"));
        let seq_length = usize::try_from(seq.length.unwrap_or_else(|| panic!("naf file had empty seq on {file:?}?"))).unwrap_or_else(|_| panic!("failed to turn u64 to usize on {file:?}"));
        results.total_length += seq_length;
        results.largest_contig = results.largest_contig.max(seq_length);
        results.shortest_contig = results.shortest_contig.min(seq_length);
        lengths.push(seq_length);
        let line = seq.sequence.unwrap_or_else(|| panic!("naf sequence had bad data {file:?}"));
        let _ = process_sequence_line(line.as_bytes(), &mut results, offset); // GC and N counts are updated
    }
    results.sequence_count = lengths.len();
    calc_nq_stats(&lengths, &mut results);

    Ok(vec![results])
}

pub fn process_gz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
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

pub fn process_zip_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
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
                    eprintln!("Error processing {file_name}: {e}");
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
            lengths.push(current_sequence_length);
            current_sequence_length = 0;
        } else {
            current_sequence_length += process_sequence_line(&line, results, offset);
        }
        line.clear();
    }

    if current_sequence_length > 0 {
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

/// Sets Everything in resuts but filename, GC count, N count
fn calc_nq_stats(lengths: &[usize], results: &mut AnalysisResults) {
    let total_length: usize = lengths.iter().sum();
    results.total_length = total_length;
    results.sequence_count = lengths.len();
    let mut cumulative_length = 0;
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));
    results.largest_contig = *sorted_lengths.first().unwrap_or(&0);
    results.shortest_contig = *sorted_lengths.last().unwrap_or(&usize::MAX);


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
