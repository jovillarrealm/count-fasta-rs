# Initial Concept
A high-performance, SIMD-accelerated tool for calculating assembly statistics from FASTA files.

---

# Product Guide: count-fasta-rs

## Overview
A high-performance, SIMD-accelerated tool for calculating assembly statistics from FASTA files. It supports various compression formats and leverages multi-threading for maximum speed and efficiency.

## Core Goals
- **Performance First:** Utilize SIMD instructions (AVX2, AVX512, NEON) and zero-copy memory mapping to maximize throughput.
- **Versatile Utility:** Provide a single, robust tool that handles a wide range of native and compressed FASTA formats without requiring external pre-processing.
- **Modern Architecture:** Leverage Rust's safety and parallel processing capabilities (`rayon`) to provide a reliable and scalable CLI tool.

## Target Audience
- **Pipeline Developers:** Focused on integrating high-performance data processing steps into automated genomic workflows.
- **Genomic Researchers:** Requiring fast and accurate assembly statistics for large-scale data analysis.

## Key Features
- **Comprehensive Format Support:** Handles plain FASTA, Gzip, XZ, Bzip2, BGZIP, Nucleotide Archive Format (NAF), and ZIP archives.
- **SIMD Acceleration:** Automatically leverages modern CPU vector extensions for nucleotide counting and sequence analysis.
- **Low Memory Footprint:** Efficiently processes massive files using `mmap` and buffered I/O.
- **Platform Ready:** Designed for easy integration into platforms like Galaxy and high-performance computing (HPC) environments.
