#!/bin/bash
set -e

# Check if cross is installed
if ! command -v cross &> /dev/null; then
    echo "Error: 'cross' is required to run ARM tests on x86."
    echo "Please install it: cargo install cross"
    exit 1
fi

echo "Running tests for aarch64-unknown-linux-gnu (NEON)..."
cross test --target aarch64-unknown-linux-gnu --verbose
