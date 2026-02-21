# count-fasta-rs

A high-performance, SIMD-accelerated tool for calculating assembly statistics from FASTA files. It supports various compression formats and leverages multi-threading for speed.

## Features

-   **Blazing Fast**: Uses SIMD instructions (AVX2, AVX512 on x86_64; NEON on AArch64) for counting nucleotides.
-   **Multi-threaded**: Processes multiple files in parallel using Rayon.
-   **Memory Efficient**: Uses memory mapping (`mmap`) and buffered reading to minimize memory footprint.
-   **Format Support**: Handles plain `.fasta`, `.fa`, `.fna` files, as well as compressed formats:
    -   Gzip (`.gz`)
    -   XZ (`.xz`)
    -   Bzip2 (`.bz2`)
    -   BGZIP (`.bgz`, `.bgzip`) via `noodles`
    -   Nucleotide Archive Format (`.naf`) via `nafcodec`
    -   ZIP archives (`.zip`) - processes compatible files inside (like the ones you get from ncbi datasets cli).

## Installation

### From Pre-built Binaries

Download the latest binary for your operating system and architecture from the [Releases page](https://github.com/jovillarrealm/count-fasta-rs/releases).

### From Source

You can build and install the latest version directly from source using Cargo:

```bash
git clone https://github.com/jovillarrealm/count-fasta-rs
cd count-fasta-rs
cargo install --path .
```

## Usage

```text
Usage: count-fasta-rs [OPTIONS] [FASTA FILE]...

Arguments:
  [FASTA FILE]...  FASTA FILE[s] to be processed [wildcards would work here].

Options:
  -c, --csv <CSV>              Path to csv to be created. It will append to the csv file if it already exists.
  -d, --directory <DIRECTORY>  Directory to be processed. Non-recursively.
  -t, --threads <THREADS>      Numbers of threads to be used. (Default: auto-detected based on CPU/files)
  -l, --legacy                 Legacy output format for debugging/compatibility.
  -h, --help                   Print help
  -V, --version                Print version
```

### Examples

**Process a single file:**
```bash
count-fasta-rs genome.fna
```

**Process multiple files using wildcards:**
```bash
count-fasta-rs *.fa.gz
```

**Process all files in a directory and save stats to a CSV:**
```bash
count-fasta-rs -d ./genomes -c stats.csv
```

**Process all files in multiple directories:**
```bash
count-fasta-rs -d ./genomes -d ./more_genomes
```

**Output format:**

Standard output:
```yaml
File name:      genome.fna 
Total length of sequence:       46759715 bp
Total number of sequences:      32
Average contig length is:       1461241 bp
Largest contig:                 18100598 bp
Shortest contig:                214 bp
N25 stats:                      25% of total sequence length is contained in the 1 sequences >= 18100598 bp
N50 stats:                      50% of total sequence length is contained in the 2 sequences >= 16068654 bp
N75 stats:                      75% of total sequence length is contained in the 3 sequences >= 10965501 bp
Total GC count:                 19534458 bp
GC %:                           41.78 %
Number of Ns:                   2900
Ns %:                           0.01 %
```

CSV output:
```csv
filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage
"genome.fna";46759715;32;1461241.09;18100598;214;16068654;41.78;2900;0.01
```

## Development & Testing

### Running Tests

To run the standard test suite:

```bash
cargo test
```

### SIMD Cross-Architecture Testing

To verify the SIMD implementations (AVX2, AVX512, NEON) across different architectures using QEMU and `cross`:

1.  **Install `cross`**: `cargo install cross`
2.  **Ensure Docker/Podman** is running.
3.  **Run the cross-test script**:

```bash
./scripts/test_simd_cross.sh
```

This script ensures correctness of the vectorized algorithms on both ARM64 and x86_64 architectures.


## License

This project is licensed under the Apache-2.0 License.
