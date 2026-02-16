---
description: Run benchmarks using hyperfine
---

Run the benchmarks using hyperfine:

1. Build the project in release mode
```bash
cargo build --release
```

// turbo
2. Run hyperfine benchmark
```bash
hyperfine -w 2 -- "target/release/count-fasta-rs -c d.csv -d bench"
```
