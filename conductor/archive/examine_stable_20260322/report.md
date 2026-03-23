# Migration Recommendation: Stable Rust with `wide` crate

## Overview
The `count-fasta-rs` tool was successfully migrated from nightly Rust (`std::simd`) to stable Rust using the `wide` crate for SIMD acceleration. This transition improves the accessibility of the tool for developers and users without significantly impacting performance.

## Key Changes
- **Rust Toolchain:** Updated `rust-toolchain.toml` to `stable`.
- **Dependencies:** Added `wide` and `bytemuck` crates.
- **SIMD Implementation:** Replaced `std::simd` with `wide::u8x32` and `wide::i8x32`.
- **Code Structure:** Refactored `src/simd.rs` to use `chunks_exact(32)` instead of the unstable `as_simd()` method.

## Benchmarking Results
Benchmarks were performed on a large genomic FASTA file (~3GB).
- **Nightly (std::simd):** 3.869 s ± 0.038 s
- **Stable (wide):** 3.947 s ± 0.095 s
- **Performance Delta:** Nightly is approximately **2% faster**, which is within the margin of error for many runs.

## Conclusion
Moving to stable Rust is highly recommended. The minimal performance loss (2%) is outweighed by the following benefits:
- **Stability:** The tool now builds on the stable channel, ensuring a more reliable build process and easier integration into CI/CD pipelines (like Galaxy).
- **Portability:** The `wide` crate provides a robust abstraction for SIMD across multiple architectures on stable Rust.
- **Ease of Use:** Users no longer need to install a nightly compiler to build the tool from source.

## Recommendation
Merge the `stable-rust-poc` branch into `main` and update the documentation to reflect the move to stable Rust.
