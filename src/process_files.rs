// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0>
// at your option. This file may not be copied, modified,
// or distributed except according to those terms.

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use liblzma::read::XzDecoder;
use memchr::memchr;
use memmap2::Mmap;
use noodles::bgzf as bgzf;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use zip::read::ZipArchive;

#[cfg(target_os = "macos")]
use std::os::unix::io::AsRawFd;
#[cfg(windows)]
use std::os::windows::fs::OpenOptionsExt;

const DECODER_BUFFER_SIZE: usize = 2 * 1024 * 1024;

pub const VALID_FILES: [&str; 3] = ["fa", "fasta", "fna"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Gzip,
    Xz,
    Bzip2,
    Bgzip,
    Zip,
    Naf,
    Unknown,
}

impl FileFormat {
    pub fn from_path(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some("gz") => FileFormat::Gzip,
            Some("xz") => FileFormat::Xz,
            Some("bz2") => FileFormat::Bzip2,
            Some("bgz") | Some("bgzip") => FileFormat::Bgzip,
            Some("zip") => FileFormat::Zip,
            Some("naf") => FileFormat::Naf,
            Some(ext) if VALID_FILES.contains(&ext) => FileFormat::Fasta,
            _ => FileFormat::Unknown,
        }
    }
}

pub fn process_any_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    match FileFormat::from_path(file) {
        FileFormat::Gzip => process_gz_file(file, no_simd, min_length),
        FileFormat::Zip => process_zip_file(file, no_simd, min_length),
        FileFormat::Xz => process_xz_file(file, no_simd, min_length),
        FileFormat::Bzip2 => process_bz2_file(file, no_simd, min_length),
        FileFormat::Bgzip => process_bgzip_file(file, no_simd, min_length),
        FileFormat::Naf => process_naf_file(file, no_simd, min_length),
        FileFormat::Fasta => process_fasta_file(file, no_simd, min_length),
        FileFormat::Unknown => Ok(Vec::new()),
    }
}

struct FastaParser {
    lengths: Vec<usize>,
    current_sequence_length: usize,
    current_gc_count: usize,
    current_n_count: usize,
    min_length: usize,
    in_header: bool,
    last_char_was_newline: bool,
    started: bool,
}

impl FastaParser {
    fn new(min_length: usize) -> Self {
        Self {
            lengths: Vec::with_capacity(250),
            current_sequence_length: 0,
            current_gc_count: 0,
            current_n_count: 0,
            min_length,
            in_header: false,
            last_char_was_newline: true, // To catch the very first '>'
            started: false,
        }
    }

    fn commit_current_contig(&mut self, results: &mut AnalysisResults) {
        if self.current_sequence_length >= self.min_length {
            if self.current_sequence_length > 0 {
                self.lengths.push(self.current_sequence_length);
            }
            results.gc_count += self.current_gc_count;
            results.n_count += self.current_n_count;
        }
        self.current_sequence_length = 0;
        self.current_gc_count = 0;
        self.current_n_count = 0;
    }

    fn feed(&mut self, data: &[u8], results: &mut AnalysisResults, no_simd: bool) {
        let mut consumed = 0;
        let len = data.len();
        while consumed < len {
            if self.in_header {
                // Find end of header
                match memchr(b'\n', &data[consumed..]) {
                    Some(pos) => {
                        consumed += pos + 1;
                        self.in_header = false;
                        self.last_char_was_newline = true;
                    }
                    None => {
                        consumed = len;
                        self.last_char_was_newline = false;
                        // Still in header for next buffer
                    }
                }
            } else {
                // In sequence
                if self.last_char_was_newline && data[consumed] == b'>' {
                    // Start of a new header
                    if self.started {
                        self.commit_current_contig(results);
                    }
                    if self.min_length == 0 {
                        results.sequence_count += 1;
                    }
                    self.started = true;
                    self.in_header = true;
                    consumed += 1;
                    self.last_char_was_newline = false;
                } else {
                    // Still in sequence, look for the next '\n>'
                    let mut next_header_start = None;
                    let mut search_pos = consumed;
                    while let Some(pos) = memchr(b'\n', &data[search_pos..]) {
                        let actual_pos = search_pos + pos;
                        if actual_pos + 1 < len {
                            if data[actual_pos + 1] == b'>' {
                                next_header_start = Some(actual_pos);
                                break;
                            }
                        } else {
                            // \n is the last char, need to check next buffer
                            break;
                        }
                        search_pos = actual_pos + 1;
                    }

                    let (chunk, chunk_end, new_last_newline) = match next_header_start {
                        Some(pos) => (&data[consumed..pos], pos + 1, true),
                        None => (&data[consumed..len], len, data[len - 1] == b'\n'),
                    };

                    if self.started {
                        let (gc, n, sc) = crate::simd::update_stats(chunk, no_simd);
                        self.current_sequence_length += sc;
                        self.current_gc_count += gc;
                        self.current_n_count += n;
                    }
                    consumed = chunk_end;
                    self.last_char_was_newline = new_last_newline;
                }
            }
        }
    }

    fn finish(mut self, results: &mut AnalysisResults) {
        if self.started {
            self.commit_current_contig(results);
        }
        if self.min_length > 0 {
            results.sequence_count = self.lengths.len();
        }
        results.calculate_stats(self.lengths);
    }
}

#[derive(Default, Clone, Debug)]
pub struct AnalysisResults {
    pub filename: String,
    pub total_length: usize,
    pub sequence_count: usize,
    pub gc_count: usize,
    pub n_count: usize,
    pub n25: usize,
    pub n25_sequence_count: usize,
    pub l25: usize,
    pub n50: usize,
    pub n50_sequence_count: usize,
    pub l50: usize,
    pub n75: usize,
    pub n75_sequence_count: usize,
    pub l75: usize,
    pub n90: usize,
    pub n90_sequence_count: usize,
    pub l90: usize,
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

    pub fn for_path(path: &Path) -> Self {
        let filename = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "unknown".to_string());
        Self::new(filename)
    }

    pub fn calculate_stats(&mut self, mut lengths: Vec<usize>) {
        if lengths.is_empty() {
            return;
        }
        let total_length: usize = lengths.iter().sum();
        self.total_length = total_length;
        self.sequence_count = lengths.len();
        let mut cumulative_length = 0;
        lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort in-place
        self.largest_contig = *lengths.first().unwrap_or(&0);
        self.shortest_contig = *lengths.last().unwrap_or(&usize::MAX);

        let threshold_25 = (total_length as u128 * 25 / 100) as usize;
        let threshold_50 = (total_length as u128 * 50 / 100) as usize;
        let threshold_75 = (total_length as u128 * 75 / 100) as usize;
        let threshold_90 = (total_length as u128 * 90 / 100) as usize;

        for (i, &length) in lengths.iter().enumerate() {
            cumulative_length += length;
            if self.n25 == 0 && cumulative_length >= threshold_25 {
                self.n25 = length;
                self.n25_sequence_count = i + 1;
                self.l25 = i + 1;
            }
            if self.n50 == 0 && cumulative_length >= threshold_50 {
                self.n50 = length;
                self.n50_sequence_count = i + 1;
                self.l50 = i + 1;
            }
            if self.n75 == 0 && cumulative_length >= threshold_75 {
                self.n75 = length;
                self.n75_sequence_count = i + 1;
                self.l75 = i + 1;
            }
            if self.n90 == 0 && cumulative_length >= threshold_90 {
                self.n90 = length;
                self.n90_sequence_count = i + 1;
                self.l90 = i + 1;
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

pub fn process_xz_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, XzDecoder::new, no_simd, min_length)
}

pub fn process_bz2_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, BzDecoder::new, no_simd, min_length)
}

pub fn process_bgzip_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, bgzf::io::Reader::new, no_simd, min_length)
}

pub fn process_fasta_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults::for_path(file);
    let file = open_file(file)?;

    // Try to use mmap, fallback to BufReader if it fails (e.g. empty file or special file)
    let mmap_result = unsafe { Mmap::map(&file) };

    match mmap_result {
        Ok(mmap) => {
            #[cfg(unix)]
            mmap.advise(memmap2::Advice::Sequential)?;

            #[cfg(target_os = "linux")]
            mmap.advise(memmap2::Advice::HugePage)?;

            process_buffer(&mmap, &mut results, no_simd, min_length)?;
        }
        Err(_) => {
            println!("Failed to mmap file: {:?}", file);
            let reader = BufReader::with_capacity(DECODER_BUFFER_SIZE, file);
            process_reader(reader, &mut results, no_simd, min_length)?;
        }
    }

    Ok(vec![results])
}

pub fn process_naf_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    let mut results = AnalysisResults::for_path(file);
    let decoder = nafcodec::Decoder::from_path(file)
        .map_err(|e| std::io::Error::other(format!("failed to open nucleotide archive: {e}")))?;

    // Process naf file
    let mut lengths = Vec::with_capacity(250);

    for may_seq in decoder {
        let seq = may_seq.map_err(|e| std::io::Error::other(format!("{file:?} had bad data: {e}")))?;
        let seq_length = usize::try_from(
            seq.length
                .ok_or_else(|| std::io::Error::other(format!("naf file had empty seq on {file:?}?")))?,
        )
        .map_err(|e| std::io::Error::other(format!("failed to turn u64 to usize on {file:?}: {e}")))?;
        if min_length > 0 && seq_length < min_length {
            continue;
        }
        lengths.push(seq_length);
        let line = seq
            .sequence
            .ok_or_else(|| std::io::Error::other(format!("naf sequence had bad data {file:?}")))?;
        let (gc, n, _) = crate::simd::update_stats(line.as_bytes(), no_simd);
        results.gc_count += gc;
        results.n_count += n;
    }
    if min_length > 0 {
        results.sequence_count = lengths.len();
    }
    results.calculate_stats(lengths);

    Ok(vec![results])
}

pub fn process_gz_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    process_decoded_stream(file, GzDecoder::new, no_simd, min_length)
}

fn process_decoded_stream<D, F>(
    file: &Path,
    decoder_factory: F,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>>
where
    D: Read,
    F: Fn(File) -> D,
{
    let mut results = AnalysisResults::for_path(file);
    let file = open_file(file)?;
    let decoder = decoder_factory(file);
    let reader = BufReader::with_capacity(DECODER_BUFFER_SIZE, decoder);
    process_reader(reader, &mut results, no_simd, min_length)?;
    Ok(vec![results])
}

pub fn process_zip_file(
    file: &Path,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<Vec<AnalysisResults>> {
    let archive_name = file
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "archive.zip".to_string());
    let file_handle = open_file(file)?;
    let buf_reader = BufReader::with_capacity(DECODER_BUFFER_SIZE, file_handle);
    let mut archive = ZipArchive::new(buf_reader)?;
    let mut all_results = Vec::new();

    for i in 0..archive.len() {
        let zip_file = archive.by_index(i)?;
        if zip_file.is_file() {
            let file_name = zip_file.name().to_owned();
            if VALID_FILES.iter().any(|&ext| file_name.ends_with(ext)) {
                let composite_name = format!("{archive_name}#{file_name}");
                let mut result = AnalysisResults::new(composite_name);
                let reader = BufReader::with_capacity(DECODER_BUFFER_SIZE, zip_file);
                if let Err(e) = process_reader(reader, &mut result, no_simd, min_length) {
                    eprintln!("Error processing {file_name}: {e}");
                    continue; // Skip this file but continue processing others
                };
                all_results.push(result);
            }
        }
    }

    Ok(all_results)
}

pub fn process_reader<R: Read>(
    mut reader: BufReader<R>,
    results: &mut AnalysisResults,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<()> {
    let mut parser = FastaParser::new(min_length);

    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        parser.feed(buf, results, no_simd);
        let consumed = buf.len();
        reader.consume(consumed);
    }

    parser.finish(results);
    Ok(())
}

#[cfg(test)]
fn update_stats(line: &[u8], results: &mut AnalysisResults, no_simd: bool) -> usize {
    let (gc, n, seq_chars) = crate::simd::update_stats(line, no_simd);
    results.gc_count += gc;
    results.n_count += n;
    seq_chars
}

pub fn process_buffer(
    data: &[u8],
    results: &mut AnalysisResults,
    no_simd: bool,
    min_length: usize,
) -> std::io::Result<()> {
    let mut parser = FastaParser::new(min_length);
    parser.feed(data, results, no_simd);
    parser.finish(results);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_calculate_stats() {
        let mut results = AnalysisResults::new("test.fa".to_string());
        let lengths = vec![100, 200, 300, 400]; // Total: 1000
        results.calculate_stats(lengths);

        assert_eq!(results.total_length, 1000);
        assert_eq!(results.sequence_count, 4);
        assert_eq!(results.largest_contig, 400);
        assert_eq!(results.shortest_contig, 100);
        assert_eq!(results.n25, 400);
        assert_eq!(results.n25_sequence_count, 1);
        assert_eq!(results.l25, 1);
        assert_eq!(results.n50, 300);
        assert_eq!(results.n50_sequence_count, 2);
        assert_eq!(results.l50, 2);
        assert_eq!(results.n75, 200);
        assert_eq!(results.n75_sequence_count, 3);
        assert_eq!(results.l75, 3);
        assert_eq!(results.n90, 200);
        assert_eq!(results.l90, 3);
        assert_eq!(results.n90_sequence_count, 3);
    }

    #[test]
    fn test_calculate_stats_single_contig() {
        let mut results = AnalysisResults::new("single.fa".to_string());
        results.calculate_stats(vec![500]);

        assert_eq!(results.total_length, 500);
        assert_eq!(results.sequence_count, 1);
        assert_eq!(results.largest_contig, 500);
        assert_eq!(results.shortest_contig, 500);
        assert_eq!(results.n25, 500);
        assert_eq!(results.l25, 1);
        assert_eq!(results.n50, 500);
        assert_eq!(results.l50, 1);
        assert_eq!(results.n75, 500);
        assert_eq!(results.l75, 1);
        assert_eq!(results.n90, 500);
        assert_eq!(results.l90, 1);
    }

    #[test]
    fn test_update_stats() {
        let mut results = AnalysisResults::default();
        update_stats(b"ATGCatgcNNnn", &mut results, false);
        assert_eq!(results.gc_count, 4);
        assert_eq!(results.n_count, 4);
    }

    #[test]
    fn test_process_buffer() {
        let data = b">seq1\nATGC\n>seq2\nAAAAA\n";
        let mut results = AnalysisResults::new("buffer".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();

        assert_eq!(results.total_length, 9);
        assert_eq!(results.sequence_count, 2);
        assert_eq!(results.gc_count, 2);
        assert_eq!(results.n_count, 0);
        assert_eq!(results.largest_contig, 5);
        assert_eq!(results.shortest_contig, 4);
    }

    #[test]
    fn test_update_stats_exhaustive() {
        let mut results = AnalysisResults::default();
        // Test all possible DNA characters in both cases, including ambiguity codes and gaps
        // Standard: AaCcGgTtNn
        // Ambiguity: RrYyWwSsMmKkHhBbVvDd
        // Gaps/Noise: - . [space]
        let input = b"AaCcGgTtNnRrYyWwSsMmKkHhBbVvDd-. \t";
        let seq_len = update_stats(input, &mut results, false);
        
        // G, g, C, c are the only 4 counted as GC
        assert_eq!(results.gc_count, 4);
        // N, n are the only 2 counted as N
        assert_eq!(results.n_count, 2);
        // Total characters excluding gaps (-.) and whitespace ( \t)
        // 10 (standard) + 20 (ambiguity) = 30
        assert_eq!(seq_len, 30);
    }

    #[test]
    fn test_process_reader_mixed_line_endings() {
        let data = b">seq1\nATGC\r\n>seq2\r\nAAAAA\n";
        let mut results = AnalysisResults::new("mixed".to_string());
        let reader = BufReader::new(&data[..]);
        process_reader(reader, &mut results, false, 0).unwrap();

        assert_eq!(results.total_length, 9);
        assert_eq!(results.sequence_count, 2);
    }

    #[test]
    fn test_process_buffer_headers_with_gc() {
        let data = b">seq_with_GC_and_N\nATGC\n>next\nNNNN\n";
        let mut results = AnalysisResults::new("headers".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();

        // Header content should NOT be counted
        assert_eq!(results.gc_count, 2); 
        assert_eq!(results.n_count, 4);
        assert_eq!(results.total_length, 8);
    }

    #[test]
    fn test_process_with_gaps_and_whitespace() {
        let data = b">seq1\nAT GC\n-..-\nATGC\n";
        let mut results = AnalysisResults::new("gaps".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();

        // ATGC (4) + ATGC (4) = 8. Gaps and spaces ignored.
        assert_eq!(results.total_length, 8);
        assert_eq!(results.gc_count, 4);
    }

    #[test]
    fn test_process_buffer_crlf() {
        let data = b">seq1\r\nATGC\r\n>seq2\r\nAAAAA\r\n";
        let mut results = AnalysisResults::new("buffer".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();

        assert_eq!(results.total_length, 9);
        assert_eq!(results.sequence_count, 2);
    }

    #[test]
    fn test_process_buffer_min_length() {
        let data = b">seq1\nAAAA\n>seq2\nAAAAAAAA\n>seq3\nAA\n"; // seq1: 4bp, seq2: 8bp, seq3: 2bp
        let mut results = AnalysisResults::new("min_len.fa".to_string());
        process_buffer(data, &mut results, false, 4).unwrap();

        assert_eq!(results.sequence_count, 2);
        assert_eq!(results.total_length, 12);
        assert_eq!(results.largest_contig, 8);
        assert_eq!(results.shortest_contig, 4);
        assert_eq!(results.n50, 8);
        assert_eq!(results.l50, 1);
    }

    #[test]
    fn test_process_buffer_min_length_all_excluded() {
        let data = b">seq1\nAAAA\n>seq2\nAA\n"; // 4bp, 2bp
        let mut results = AnalysisResults::new("min_len_all.fa".to_string());
        process_buffer(data, &mut results, false, 10).unwrap();

        assert_eq!(results.sequence_count, 0);
        assert_eq!(results.total_length, 0);
        assert_eq!(results.largest_contig, 0);
        assert_eq!(results.shortest_contig, usize::MAX);
        assert_eq!(results.n50, 0);
        assert_eq!(results.l50, 0);
    }

    #[test]
    fn test_process_empty() {
        let data = b"";
        let mut results = AnalysisResults::new("empty".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();
        assert_eq!(results.total_length, 0);
        assert_eq!(results.sequence_count, 0);

        let mut results2 = AnalysisResults::new("empty_reader".to_string());
        let reader = BufReader::new(&data[..]);
        process_reader(reader, &mut results2, false, 0).unwrap();
        assert_eq!(results2.total_length, 0);
        assert_eq!(results2.sequence_count, 0);
    }

    #[test]
    fn test_process_only_header() {
        let data = b">only_header\n";
        let mut results = AnalysisResults::new("only_header".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();
        assert_eq!(results.total_length, 0);
        assert_eq!(results.sequence_count, 1);

        let mut results2 = AnalysisResults::new("only_header_reader".to_string());
        let reader = BufReader::new(&data[..]);
        process_reader(reader, &mut results2, false, 0).unwrap();
        assert_eq!(results2.total_length, 0);
        assert_eq!(results2.sequence_count, 1);
    }

    #[test]
    fn test_process_no_trailing_newline() {
        let data = b">seq1\nATGC";
        let mut results = AnalysisResults::new("no_newline".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();
        assert_eq!(results.total_length, 4);
        assert_eq!(results.sequence_count, 1);

        let mut results2 = AnalysisResults::new("no_newline_reader".to_string());
        let reader = BufReader::new(&data[..]);
        process_reader(reader, &mut results2, false, 0).unwrap();
        assert_eq!(results2.total_length, 4);
        assert_eq!(results2.sequence_count, 1);
    }

    #[test]
    fn test_process_lines_before_header() {
        let data = b"some noise\n>seq1\nATGC\n";
        let mut results = AnalysisResults::new("noise".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();
        // Noise is now correctly ignored.
        assert_eq!(results.sequence_count, 1);
        assert_eq!(results.total_length, 4);
    }

    #[test]
    fn test_real_world_complexities() {
        let data = b"; legacy comment line\n>seq1 with spaces\nATGC\n>seq1\nAAAA\n>  seq2\tmetadata\nGGGG\n";
        let mut results = AnalysisResults::new("complex".to_string());
        process_buffer(data, &mut results, false, 0).unwrap();

        // 1. Comment line is ignored. 3 sequences found.
        // 2. Total length: 4 (ATGC) + 4 (AAAA) + 4 (GGGG) = 12
        assert_eq!(results.total_length, 12);
        assert_eq!(results.sequence_count, 3); 
    }

    #[test]
    fn test_process_missing_file() {
        let path = Path::new("non_existent_file.fa");
        let res = process_fasta_file(path, false, 0);
        assert!(res.is_err());
    }

    #[test]
    fn test_process_reader_io_error() {
        struct FailingReader;
        impl Read for FailingReader {
            fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
                Err(std::io::Error::new(std::io::ErrorKind::Other, "forced failure"))
            }
        }
        let reader = BufReader::new(FailingReader);
        let mut results = AnalysisResults::new("failing".to_string());
        let res = process_reader(reader, &mut results, false, 0);
        assert!(res.is_err());
    }

    #[test]
    fn test_process_corrupted_gz() {
        let mut temp_file = std::env::temp_dir();
        temp_file.push("corrupted.gz");
        fs::write(&temp_file, b"this is not a valid gzip file").unwrap();
        let res = process_gz_file(&temp_file, false, 0);
        assert!(res.is_err());
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_process_corrupted_xz() {
        let mut temp_file = std::env::temp_dir();
        temp_file.push("corrupted.xz");
        fs::write(&temp_file, b"this is not a valid xz file").unwrap();
        let res = process_xz_file(&temp_file, false, 0);
        assert!(res.is_err());
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_process_corrupted_zip() {
        let mut temp_file = std::env::temp_dir();
        temp_file.push("corrupted.zip");
        fs::write(&temp_file, b"this is not a valid zip file").unwrap();
        let res = process_zip_file(&temp_file, false, 0);
        assert!(res.is_err());
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_mmap_fallback_empty_file() {
        let mut temp_file = std::env::temp_dir();
        temp_file.push("empty_mmap.fa");
        fs::write(&temp_file, b"").unwrap();
        let res = process_fasta_file(&temp_file, false, 0);
        assert!(res.is_ok());
        let results = res.unwrap();
        assert_eq!(results[0].sequence_count, 0);
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_process_invalid_naf() {
        let mut temp_file = std::env::temp_dir();
        temp_file.push("invalid.naf");
        fs::write(&temp_file, b"this is not a valid naf file").unwrap();
        let res = process_naf_file(&temp_file, false, 0);
        assert!(res.is_err());
        let _ = fs::remove_file(temp_file);
    }
}
