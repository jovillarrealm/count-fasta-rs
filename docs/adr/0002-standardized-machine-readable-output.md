# 2. Standardized Machine-Readable Output Formats

## Context

The legacy implementation emitted semicolon-separated values under `-c/--csv` without support for standard RFC 4180 comma delimiters or alternative serialization formats such as TSV and JSON. Furthermore, reports omitted critical metrics such as N90, L50, and L90 from tabular exports. Modern bioinformatics pipelines (Nextflow, Snakemake, Python Pandas, R) require standardized delimiters and structured formats for automated ingestion.

## Decision

To avoid breaking existing users' workflows and legacy scripts, we preserve backwards compatibility by default:
- Default `-c, --csv <FILE>` continues to emit semicolon-delimited values matching historical output and schemas.
- Default terminal stdout output continues to display the classic metrics suite (including N50).
- Standard RFC 4180 comma-separated values (CSV) can be requested using the explicit `--comma` flag.
- Extended continuity metrics ($N_{25}, N_{75}, N_{90}$ and $L_{25}, L_{50}, L_{75}, L_{90}$) can be requested across interactive and machine-readable exports using the `-e, --extended` flag.
- Dedicated `--tsv` and `--json` format options provide native TSV and JSON serialization.

## Consequences

- Existing user scripts, Excel workflows, and regression tests relying on semicolon-delimited CSV and classic terminal output continue to run without any breakage.
- Modern bioinformatics pipelines (Nextflow, Snakemake, Pandas, R) can opt into RFC 4180 comma-separated records via `--comma`.
- Detailed continuity metrics ($N_{25}$–$N_{90}$ and $L_{25}$–$L_{90}$) are available on demand via `-e, --extended` across terminal, CSV, TSV, and JSON formats.
- The `--legacy` flag is retained for explicit historical debugging.
