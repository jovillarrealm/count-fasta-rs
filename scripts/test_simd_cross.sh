#!/bin/bash
set -e

# Ensure cross is installed
if ! command -v cross &> /dev/null; then
    echo "Error: 'cross' is not installed. Run: cargo install cross"
    exit 1
fi

echo "Running SIMD consistency tests for multiple architectures..."

# 1. Test AArch64 (NEON)
echo "----------------------------------------------------"
echo "Testing AArch64 (NEON) via QEMU..."
# QEMU_CPU=max ensures NEON and other features are enabled in the emulator
# We use simd::tests to run the specific SIMD consistency tests
QEMU_CPU=max cross test --target aarch64-unknown-linux-gnu --bin count-fasta-rs simd::tests

# 2. Test x86_64 (AVX2/AVX512)
echo "----------------------------------------------------"
echo "Testing x86_64 (AVX2/AVX512)..."
# If host is x86_64, this runs natively. 
# If host is AArch64 (e.g. Apple Silicon), this runs via QEMU.
# Note: To test AVX512 on a non-AVX512 x86_64 host, you might need a custom runner 
# or run on an AArch64 host via QEMU with -cpu max.
QEMU_CPU=max cross test --target x86_64-unknown-linux-gnu --bin count-fasta-rs simd::tests

echo "----------------------------------------------------"
echo "All SIMD tests completed successfully!"
