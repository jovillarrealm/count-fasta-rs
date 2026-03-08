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

You can build and install the latest version directly from source using Cargo. Note that this project requires the **Nightly** Rust channel for portable SIMD support:

```bash
git clone https://github.com/jovillarrealm/count-fasta-rs
cd count-fasta-rs
cargo +nightly install --path .
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

## Architecture & Performance

`count-fasta-rs` achieves its performance through three main architectural pillars:

- **Zero-Copy Memory Mapping**: Plain text FASTA files are memory-mapped (`mmap`) into the application's address space. The application uses OS-level hints (like `madvise` with `MADV_SEQUENTIAL` and `MADV_HUGEPAGE`) to bypass userspace buffering entirely. This allows ingestion speeds to scale up to the physical limits of the underlying storage.
- **Portable SIMD**: The sequence parsing and nucleotide counting logic leverages Rust's Nightly `#![feature(portable_simd)]`. This provides a single, safe, and highly maintainable codebase that automatically compiles down to heavily optimized vector instructions (e.g., AVX2, AVX-512, NEON, SVE2) depending on the target hardware without relying on brittle `core::arch` intrinsics.
- **Multiprocessing**: The application processes files concurrently using the `rayon` crate, allocating one file per CPU core for linear scaling across massive datasets.

## Development & Testing Stack

To fully develop, test, and release this application, the following software stack is required:

### 1. Standard Testing
Run the standard native test suite:
```bash
cargo test
```

### 2. Cross-Architecture Testing (QEMU)
To test SIMD implementations on non-native architectures (e.g., testing ARM NEON from an x86_64 Linux host), the project uses QEMU:
- **`cross`**: Rust cross-compilation tool (`cargo install cross`).
- **Docker or Podman**: Container runtime for `cross` to execute within.
- **`qemu-user-static`**: Must be installed on the host machine to allow Docker to run foreign architectures. The test script explicitly passes `QEMU_CPU=max` to ensure all vector extensions are emulated.

Run the cross-tests:
```bash
./scripts/test_simd_cross.sh
```

### 3. Performance Benchmarking
To verify performance regressions, a dedicated script compares the Rust implementation against the legacy Perl script:
- **`hyperfine`**: Command-line benchmarking tool (`cargo install hyperfine`).
- **`perl`**: Required to run the legacy baseline (`alternatives/count_fasta_cnsg.pl`).

Run the benchmarks (requires files in the `bench/` directory):
```bash
./scripts/benchmark.sh
```
Results are saved as Markdown reports in the `performance/` directory.

### 4. Releasing (GoReleaser)
Release automation is handled by GoReleaser, building cross-platform binaries natively or via Zig:
- **GoReleaser**: The release automation tool.
- **`cross`**: For building Linux and Windows targets in containers.
- **`cargo-zigbuild`**: For compiling macOS (Apple Darwin) targets cleanly from Linux (`cargo install --locked cargo-zigbuild`).
- **Docker/Podman**: Required by `cross` during the release process.

## License

This project is licensed under the Apache-2.0 License.
