# Tech Stack: count-fasta-rs

## Core Language
- **Rust (Nightly, Edition 2024):** Required for `#![feature(portable_simd)]` and modern language features.

## Key Libraries & Frameworks
- **CLI Handling:** `clap` (v4 with derive feature) for robust command-line argument parsing.
- **Parallelism:** `rayon` for data-parallel processing of multiple files.
- **Compression Support:**
    - `flate2` (Gzip)
    - `bzip2` (Bzip2)
    - `liblzma` (XZ)
    - `noodles` (BGZIP support)
    - `nafcodec` (Nucleotide Archive Format support)
    - `zip` (Archive file processing)
- **Performance Optimizations:**
    - `memmap2` (Zero-copy memory mapping for FASTA files)
    - `memchr` (Highly optimized byte search)
    - `num_cpus` (Automatic thread pool scaling)
- **SIMD Support:** `std::simd` (via nightly portable-simd feature) for hardware-accelerated nucleotide counting.

## Build & CI/CD
- **Testing:** `cargo test` for unit and integration tests.
- **Cross-Arch Testing:** `cross` and QEMU for ARM/x86_64 verification.
- **Performance Benchmarking:** `hyperfine` for regression testing against Perl baselines.
- **Release Automation:** `GoReleaser` with `cargo-zigbuild` for cross-platform binaries.
