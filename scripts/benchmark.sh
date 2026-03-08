#!/bin/bash
set -e

# Ensure we have hyperfine installed
if ! command -v hyperfine &> /dev/null; then
    echo "Error: hyperfine is not installed. Please install it first."
    exit 1
fi

echo "Building release binary..."
cargo build --release

BIN="./target/release/count-fasta-rs"
PERL_SCRIPT="./alternatives/count_fasta_cnsg.pl"

if [ ! -f "$PERL_SCRIPT" ]; then
    echo "Error: Perl script not found at $PERL_SCRIPT"
    exit 1
fi

# Find the largest .fa or .fna file in bench/
# Using stat to get size and filename (%n), sort numerically, pick largest, extract filename
LARGEST_FILE=$(find bench/ -maxdepth 1 -type f \( -name "*.f*a" \) -exec stat -c "%s %n" {} + | sort -n | tail -1 | awk '{print $2}')

if [ -z "$LARGEST_FILE" ]; then
    echo "Error: No .fa or .fna files found in bench/"
    exit 1
fi

mkdir -p performance

echo "=========================================================="
echo "Benchmarking Largest File: $LARGEST_FILE"
echo "=========================================================="

hyperfine --warmup 3 --export-markdown performance/largest_file_results.md --show-output \
    -n "Perl" "perl $PERL_SCRIPT '$LARGEST_FILE'" \
    -n "Rust (No SIMD)" "$BIN --no-simd '$LARGEST_FILE'" \
    -n "Rust (SIMD)" "$BIN '$LARGEST_FILE'"

echo "=========================================================="
echo "Benchmarking Entire Directory: bench/"
echo "=========================================================="

# Find all valid files to pass to the Perl script, since it takes multiple arguments
# Joining files with spaces for the command line
ALL_FILES=$(find bench/ -maxdepth 1 -type f \( -name "*.f*a" \) | paste -sd " " -)

hyperfine --warmup 1 --export-markdown performance/all_files_results.md \
    -n "Perl" "for f in $ALL_FILES; do perl $PERL_SCRIPT \"\$f\"; done" \
    -n "Rust (No SIMD, Single-threaded)" "$BIN --no-simd -t 1 -d bench/" \
    -n "Rust (No SIMD, Multi-threaded)" "$BIN --no-simd -d bench/" \
    -n "Rust (SIMD, Single-threaded)" "$BIN -t 1 -d bench/" \
    -n "Rust (SIMD, Multi-threaded)" "$BIN -d bench/"

echo "Benchmarks completed. Results saved to performance/ directory."
