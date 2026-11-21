# count-fasta-rs Domain Glossary

Core domain model and canonical terminology for `count-fasta-rs`, a high-performance assembly statistics calculator.

## Language

**Assembly**:
A reconstructed set of genomic sequences representing an organism's or community's genome, evaluated as a whole for fragmentation, size, and composition.
_Avoid_: Read set, alignment, sequencing run

**Contig**:
A contiguous reconstructed sequence within an assembly, corresponding to a single FASTA sequence entry in this tool.
_Avoid_: Read, transcript, k-mer

**FASTA Record**:
A discrete header line (prefixed by `>`) followed by one or more lines of nucleotide sequence characters.
_Avoid_: FASTA blob, entry chunk

**Assembly Metrics**:
A suite of statistical summarizations characterizing assembly continuity, size distribution, and nucleotide composition (total length, contig counts, average length, largest/shortest contigs, N50, GC%, and N% by default, with extended continuity metrics N25/N75/N90 and L25/L50/L75/L90 available via `-e, --extended`).
_Avoid_: Quality scores, variant statistics, alignment summary

**Nx Metric (N25, N50, N75, N90)**:
The contig length \(L\) such that sequences of length \(\ge L\) represent at least \(x\%\) of the total assembly length.
_Avoid_: Median contig length, average length

**Lx Metric (L25, L50, L75, L90)**:
The minimum count of contigs whose combined length represents at least \(x\%\) of the total assembly length.
_Avoid_: Contig rank, sequence count threshold

**Length Threshold**:
A user-specified minimum contig length cutoff (\(\ge L\)) filtering which contigs are included in the calculated assembly metrics.
_Avoid_: Size filter, contig trimmer

**Assembly-Level Aggregation**:
Summarization performed strictly across the whole assembly, producing one aggregate result per input assembly rather than per-contig breakdowns.
_Avoid_: Per-contig reporting, contig dump

**Archive Member Identifier**:
A composite identifier string (`archive_file#internal_path`) uniquely establishing provenance for sequences analyzed from within compressed archives (e.g. ZIP).
_Avoid_: Plain basename, zip entry

**Piped Stream**:
A sequential data stream supplied via standard input (`stdin`) rather than a direct filesystem path or archive entry.
_Avoid_: Input descriptor, raw stream

**Resilient Batch Execution**:
A fault-tolerant processing strategy where single file read or parse errors are logged to `stderr` while valid files continue to be processed, culminating in a non-zero exit code if any error occurred.
_Avoid_: Fail-fast, silent fail
