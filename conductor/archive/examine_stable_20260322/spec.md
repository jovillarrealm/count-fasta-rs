# Specification: Examine if this can be moved to latest stable

## Overview
Currently, the `count-fasta-rs` tool relies on the Rust nightly toolchain specifically for the `portable_simd` feature. This track aims to evaluate if the tool can be migrated to the latest stable Rust version (Edition 2024 or later) to improve stability and ease of installation for end-users, especially in environments like Galaxy.

## Goals
- Determine the current status of `std::simd` in stable Rust.
- Identify and test alternative SIMD crates (e.g., `wide`, `packed_simd`, or architecture-specific intrinsics) if `std::simd` remains unstable.
- Assess the performance impact of moving away from `portable_simd`.
- Propose a migration path to stable Rust if feasible.

## Acceptance Criteria
- A comprehensive report on the feasibility of moving to stable Rust.
- A proof-of-concept (PoC) implementation on a branch demonstrating the tool running on stable Rust (if feasible).
- Performance benchmarks comparing the nightly and stable versions.

## Out of Scope
- Rewriting the entire core logic unless necessary for the stable migration.
- Implementing support for new FASTA formats during this track.
