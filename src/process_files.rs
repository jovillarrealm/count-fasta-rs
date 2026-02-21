// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use liblzma::read::XzDecoder;
use memchr::memchr;
use memmap2::Mmap;
use noodles_bgzf as bgzf;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use zip::read::ZipArchive;

#[cfg(target_os = "macos")]
use std::os::unix::io::AsRawFd;
#[cfg(windows)]
use std::os::windows::fs::OpenOptionsExt;

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

impl AnalysisResults {
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            shortest_contig: usize::MAX,
            ..Default::default()
        }
    }

    pub fn calculate_stats(&mut self, mut lengths: Vec<usize>) {
        let total_length: usize = lengths.iter().sum();
        self.total_length = total_length;
        self.sequence_count = lengths.len();
        let mut cumulative_length = 0;
        lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort in-place
        self.largest_contig = *lengths.first().unwrap_or(&0);
        self.shortest_contig = *lengths.last().unwrap_or(&usize::MAX);

        for (i, &length) in lengths.iter().enumerate() {
            cumulative_length += length;
            if self.n25 == 0 && cumulative_length >= total_length / 4 {
                self.n25 = length;
                self.n25_sequence_count = i + 1;
            }
            if self.n50 == 0 && cumulative_length >= total_length / 2 {
                self.n50 = length;
                self.n50_sequence_count = i + 1;
            }
            if self.n75 == 0 && cumulative_length >= total_length * 3 / 4 {
                self.n75 = length;
                self.n75_sequence_count = i + 1;
                break;
            }
        }
    }
}

pub fn open_file<P: AsRef<Path>>(path: P) -> std::io::Result<File> {
    let mut options = OpenOptions::new();
    options.read(true);

    #[cfg(windows)]
    {
        const FILE_FLAG_SEQUENTIAL_SCAN: u32 = 0x08000000;
        options.attributes(FILE_FLAG_SEQUENTIAL_SCAN);
    }

    let file = options.open(path)?;

    #[cfg(target_os = "macos")]
    {
        unsafe extern "C" {
            fn fcntl(
                fd: std::os::raw::c_int,
                cmd: std::os::raw::c_int,
                arg: std::os::raw::c_int,
            ) -> std::os::raw::c_int;
        }
        const F_RDAHEAD: std::os::raw::c_int = 45;
        unsafe {
            fcntl(file.as_raw_fd(), F_RDAHEAD, 1);
        }
    }

    Ok(file)
}

pub fn process_xz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, buffer_size, XzDecoder::new)
}

pub fn process_bz2_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, buffer_size, BzDecoder::new)
}

pub fn process_bgzip_file(
    file: &Path,
    buffer_size: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, buffer_size, bgzf::io::Reader::new)
}

pub fn process_fasta_file(
    file: &Path,
    buffer_size: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults::new(file.file_name().unwrap().to_string_lossy().to_string());
    let file = open_file(file)?;

    // Try to use mmap, fallback to BufReader if it fails (e.g. empty file or special file)
    let mmap_result = unsafe { Mmap::map(&file) };

    match mmap_result {
        Ok(mmap) => {
            #[cfg(unix)]
            mmap.advise(memmap2::Advice::Sequential)?;

            #[cfg(target_os = "linux")]
            mmap.advise(memmap2::Advice::HugePage)?;

            process_buffer(&mmap, &mut results)?;
        }
        Err(_) => {
            // Re-open or just use the existing file handle if possible,
            // but File definition in map consumes it? No, map takes &File.
            // But we need to reset position if we read from it?
            // We haven't read from it yet.
            println!("Failed to mmap file: {:?}", file);
            let reader = BufReader::with_capacity(buffer_size, file);
            process_reader(reader, &mut results)?;
        }
    }

    Ok(vec![results])
}

pub fn process_naf_file(file: &Path) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults::new(file.file_name().unwrap().to_string_lossy().to_string());
    let decoder = nafcodec::Decoder::from_path(file).expect("failed to open nucleotide archive");

    // Process naf file
    let mut lengths = Vec::with_capacity(250);

    for may_seq in decoder {
        let seq = may_seq.unwrap_or_else(|_| panic!("{file:?} had bad data"));
        let seq_length = usize::try_from(
            seq.length
                .unwrap_or_else(|| panic!("naf file had empty seq on {file:?}?")),
        )
        .unwrap_or_else(|_| panic!("failed to turn u64 to usize on {file:?}"));
        results.total_length += seq_length;
        results.largest_contig = results.largest_contig.max(seq_length);
        results.shortest_contig = results.shortest_contig.min(seq_length);
        lengths.push(seq_length);
        let line = seq
            .sequence
            .unwrap_or_else(|| panic!("naf sequence had bad data {file:?}"));
        update_stats(line.as_bytes(), &mut results);
    }
    results.sequence_count = lengths.len();
    results.calculate_stats(lengths);

    Ok(vec![results])
}

pub fn process_gz_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, buffer_size, GzDecoder::new)
}

fn process_decoded_stream<D, F>(
    file: &Path,
    buffer_size: usize,
    decoder_factory: F,
) -> std::io::Result<Vec<AnalysisResults>>
where
    D: Read,
    F: Fn(File) -> D,
{
    let mut results = AnalysisResults::new(file.file_name().unwrap().to_string_lossy().to_string());
    let file = open_file(file)?;
    let decoder = decoder_factory(file);
    let reader = BufReader::with_capacity(buffer_size, decoder);
    process_reader(reader, &mut results)?;
    Ok(vec![results])
}

pub fn process_zip_file(file: &Path, buffer_size: usize) -> std::io::Result<Vec<AnalysisResults>> {
    let file = open_file(file)?;
    let buf_reader = BufReader::with_capacity(buffer_size, file);
    let mut archive = ZipArchive::new(buf_reader)?;
    let mut all_results = Vec::new();

    for i in 0..archive.len() {
        let zip_file = archive.by_index(i)?;
        if zip_file.is_file() {
            let file_name = zip_file.name().to_owned();
            if VALID_FILES.iter().any(|&ext| file_name.ends_with(ext)) {
                let mut result = AnalysisResults::new(
                    Path::new(&file_name)
                        .file_name()
                        .unwrap()
                        .to_string_lossy()
                        .to_string(),
                );
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

    // Read header line safely to determine offset
    let mut line = Vec::new(); // Temporary buffer for header only
    let offset = if reader.read_until(b'\n', &mut line)? > 0 {
        results.sequence_count += 1;
        if line.ends_with(b"\r\n") {
            Some(2)
        } else if line.ends_with(b"\n") {
            Some(1)
        } else {
            None
        }
    } else {
        return Ok(());
    };

    let Some(offset) = offset else {
        return Ok(());
    };

    let mut start_of_line = true;
    let mut in_header = false;

    loop {
        let buf = reader.fill_buf()?;
        let len = buf.len();
        if len == 0 {
            break;
        }

        let mut consumed = 0;

        while consumed < len {
            // Check start of header
            if start_of_line && buf[consumed] == b'>' {
                lengths.push(current_sequence_length);
                current_sequence_length = 0;
                in_header = true;
            }

            if in_header {
                // Scan for end of header (newline)
                match memchr(b'\n', &buf[consumed..]) {
                    Some(pos) => {
                        consumed += pos + 1;
                        start_of_line = true;
                        in_header = false;
                    }
                    None => {
                        consumed = len;
                        start_of_line = false; // Still in header
                    }
                }
            } else {
                // In sequence
                match memchr(b'\n', &buf[consumed..]) {
                    Some(pos) => {
                        let end = consumed + pos + 1;
                        let chunk = &buf[consumed..end];
                        update_stats(chunk, results);
                        current_sequence_length += chunk.len().saturating_sub(offset);

                        consumed = end;
                        start_of_line = true;
                    }
                    None => {
                        // Entire buffer is sequence data (no newline)
                        let chunk = &buf[consumed..];
                        update_stats(chunk, results);
                        current_sequence_length += chunk.len();

                        consumed = len;
                        start_of_line = false;
                    }
                }
            }
        }
        reader.consume(consumed);
    }

    if current_sequence_length > 0 {
        lengths.push(current_sequence_length);
    }

    results.calculate_stats(lengths);
    Ok(())
}

fn update_stats(line: &[u8], results: &mut AnalysisResults) {
    results.gc_count += bytecount::count(line, b'G')
        + bytecount::count(line, b'g')
        + bytecount::count(line, b'C')
        + bytecount::count(line, b'c');
    results.n_count += bytecount::count(line, b'N') + bytecount::count(line, b'n');
}

fn process_buffer(data: &[u8], results: &mut AnalysisResults) -> std::io::Result<()> {
    let mut lengths = Vec::with_capacity(250);
    let mut current_sequence_length = 0;

    let mut i = 0;
    let len = data.len();

    // Check if empty
    if len == 0 {
        return Ok(());
    }

    // find first newline
    let first_newline = match memchr(b'\n', &data[i..]) {
        Some(pos) => pos + i,
        None => return Ok(()), // No newline found
    };

    let first_line = &data[i..=first_newline];
    results.sequence_count += 1;

    let offset = if first_line.ends_with(b"\r\n") {
        2
    } else if first_line.ends_with(b"\n") {
        1
    } else {
        0
    };

    i = first_newline + 1;

    while i < len {
        let next_newline = match memchr(b'\n', &data[i..]) {
            Some(pos) => pos + i,
            None => len - 1,
        };

        let line_end = if next_newline < len {
            next_newline + 1
        } else {
            len
        };
        let line = &data[i..line_end];

        if line.first() == Some(&b'>') {
            lengths.push(current_sequence_length);
            current_sequence_length = 0;
            // Start new sequence
        } else {
            update_stats(line, results);
            let len = line.len();
            current_sequence_length += if line.ends_with(b"\n") {
                len - offset
            } else {
                len
            };
        }

        i = line_end;
    }

    if current_sequence_length > 0 {
        lengths.push(current_sequence_length);
    }

    results.calculate_stats(lengths);
    Ok(())
}
