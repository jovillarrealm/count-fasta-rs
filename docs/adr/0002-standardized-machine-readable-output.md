# 2. Standardized Machine-Readable Output Formats

## Context

The legacy implementation emitted semicolon-separated values under `-c/--csv` without support for standard RFC 4180 comma delimiters or alternative serialization formats such as TSV and JSON. Furthermore, reports omitted critical metrics such as N90, L50, and L90 from tabular exports. Modern bioinformatics pipelines (Nextflow, Snakemake, Python Pandas, R) require standardized delimiters and structured formats for automated ingestion.

## Decision

We standardize machine-readable exports to RFC 4180-compliant comma-separated values (CSV) for `--csv`, while offering explicit `--tsv` and `--json` format options. All metrics—including N25, N50, N75, N90, L25, L50, L75, and L90—are unified and output by default without requiring an extra configuration flag.

## Consequences

- Direct interoperability with standard CSV/TSV parsers without custom delimiter flags.
- Comprehensive assembly continuity statistics (N25–N90 and L25–L90) are always available in both interactive and tabular/JSON outputs without multi-tiered flags.
- Legacy semicolon output is preserved only under the `--legacy` debugging flag.
