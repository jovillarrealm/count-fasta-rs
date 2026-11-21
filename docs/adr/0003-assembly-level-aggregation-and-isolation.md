# 3. Assembly-Level Aggregation and File Isolation

## Context

Users frequently ask bioinformatics tools to either dump per-contig metric tables (millions of rows for fragmented assemblies) or merge multi-file inputs (such as separate chromosome files) into a single composite assembly statistic. We evaluated whether to support per-contig breakdowns and cross-file aggregation.

## Decision

`count-fasta-rs` operates strictly at the whole-assembly level and evaluates each file as an isolated, independent assembly. It will not generate per-contig output tables, nor will it merge multiple files into a single cross-file assembly summary.

## Consequences

- Rayon parallelism across input files remains completely lock-free and embarrassingly parallel without inter-thread synchronization.
- Output serialization remains lightweight, producing exactly one aggregate metric record per input file or archive member.
- Workflows needing per-contig data are directed to specialized per-record tools (`samtools faidx`, `seqkit fx2tab`).
