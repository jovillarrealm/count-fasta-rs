#!/bin/bash

# Check if cargo-llvm-cov is installed
if ! command -v cargo-llvm-cov &> /dev/null; then
    echo "cargo-llvm-cov is not installed. You can install it with:"
    echo "cargo install cargo-llvm-cov"
    exit 1
fi

echo "Running code coverage..."
cargo llvm-cov --all-features --workspace --lcov --output-path lcov.info
cargo llvm-cov report --html
echo "Coverage report generated in target/llvm-cov/html/index.html"
