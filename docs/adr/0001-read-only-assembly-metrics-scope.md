# 1. Strictly Read-Only Assembly Metrics Scope

## Context

Many bioinformatics sequence processing tools (e.g., `seqkit`, `fastp`, `quast`) expand beyond metrics into sequence manipulation, filtering, adapter trimming, format conversion, or reference-based alignment. We needed to define the fundamental scope boundaries of `count-fasta-rs`.

## Decision

`count-fasta-rs` is strictly a read-only statistics and assembly assessment engine. It will never modify sequences, filter records, or perform reference-based alignments.

## Consequences

- Vectorized SIMD parsing routines can focus exclusively on counting nucleotides and delimiters without supporting buffering for transformations or rewrites.
- Memory consumption remains minimal across all supported formats.
- Downstream tasks like sequence filtering or polishing are explicitly delegated to specialized pipeline tools.
