---
description: Profile the application using samply
---

To identify bottlenecks in `count-fasta-rs` using `samply`:

1.  **Build with profiling support:**
    The project has a `profiling` profile that inherits from `release` but keeps debug symbols.
    ```bash
    cargo build --profile profiling
    ```

2.  **Record a profiling session:**
    Run `samply record` followed by the command you want to profile. For example, processing the human genome files in the `bench` directory:
    ```bash
    samply record target/profiling/count-fasta-rs -d bench
    ```
    *Note: If you are processing many files, you might want to limit it to one or two to keep the profile trace manageable.*

3.  **Analyze the results:**
    `samply` will automatically open a web-based interface (using Firefox Profiler) once the command finishes.
    - Look at the **Call Tree** to see where most of the time is spent.
    - Use the **Flame Graph** to visualize the stack traces.
    - Look for `process_reader`, `process_buffer`, or `update_stats` to see their relative costs.

4.  **Common Bottlenecks to look for:**
    - **I/O Wait:** If the CPU is idle while waiting for disk.
    - **`update_stats`:** If the character counting logic is taking a large percentage of time.
    - **`memchr`:** Overhead of finding newlines.
    - **Decompression:** If processing `.gz` or `.bz2` files, the decoder might be the bottleneck.
