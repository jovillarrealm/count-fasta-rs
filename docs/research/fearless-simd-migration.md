# Research: Migrating from `wide` to `fearless_simd`

## Executive Summary

`count-fasta-rs` currently relies on `wide = "1.2.0"` and `bytemuck = "1.25.0"` in `src/simd.rs` using fixed 32-byte chunks (`u8x32`). While functional, this restricts processing to 256-bit vectors, requiring emulation on 128-bit architectures (e.g. ARM NEON) and leaving 512-bit vector registers underutilized on AVX-512 hardware.

In September 2026, the Linebender organization released **`fearless_simd` v1.0.0** (MSRV 1.89). `fearless_simd` provides safe, zero-dependency, portable SIMD with runtime dynamic multiversioning (`dispatch!`) across SSE2, AVX2, AVX-512, and ARM NEON without requiring `unsafe` or nightly compiler features.

Migrating `count-fasta-rs` to `fearless_simd` will:
1. Provide true hardware-native SIMD vector widths (16-byte on NEON/SSE, 32-byte on AVX2, 64-byte on AVX-512).
2. Eliminate `bytemuck` casts for bitmask extraction (native `mask.to_bitmask().count_ones()`).
3. Maintain 100% safe Rust without nightly `portable_simd` compiler flags.

---

## Primary Source Findings

### 1. Crate Overview & Dependencies
- **Crate**: `fearless_simd` v1.0.0 (crates.io: [fearless_simd](https://crates.io/crates/fearless_simd/1.0.0), docs: [docs.rs/fearless_simd](https://docs.rs/fearless_simd/1.0.0)).
- **Companion Macro Crate**: `fearless_simd_macros = "0.1.0"` provides the `#[simd]` attribute for automatic multi-architecture function generation.
- **MSRV**: Rust 1.89 (Current project compiler is `rustc 1.99.0`).
- **Dependencies**: The core `fearless_simd` crate has zero external dependencies (`default = ["std"]`).

### 2. Architecture & Vector Abstraction
- **Dynamic Multiversioning (`dispatch!`)**:
  `Level::new()` detects CPU features at runtime. The `dispatch!(level, simd => function(simd, ...))` macro automatically selects and executes the optimal implementation (e.g., AVX-512, AVX2, SSE2, or NEON) without manual CPU feature queries.
- **Hardware-Native Widths**:
  Instead of hardcoding 32-byte vectors (`u8x32`), vector types are accessed via the `Simd` capability trait:
  ```rust
  let chunk_len = S::u8s::LEN; // 16 on SSE/NEON, 32 on AVX2, 64 on AVX-512
  ```
- **Vector Operations**:
  - `S::u8s::from_slice(simd, chunk)`: Loads a slice into the native SIMD vector.
  - `v.simd_eq(target)`: Produces a vector mask.
  - `mask.to_bitmask()`: Returns a scalar bitmask (integer) representing lane booleans.
  - `.count_ones()`: Directly counts matched lanes without `bytemuck::cast` workarounds.

---

## Comparison: `wide` vs `fearless_simd`

| Dimension | Current (`wide` 1.2.0) | Proposed (`fearless_simd` 1.0.0) |
| :--- | :--- | :--- |
| **Vector Width** | Fixed 32 bytes (`u8x32`) | Dynamic native width (`S::u8s::LEN`: 16, 32, 64 bytes) |
| **AVX-512 Support** | No (limited to 256-bit AVX2) | Yes (full 512-bit vectorization) |
| **AArch64 / NEON** | Emulated / split across 128-bit | Native 128-bit NEON execution |
| **Bitmask Counting** | Requires `bytemuck::cast` to `i8x32` | Built-in `mask.to_bitmask().count_ones()` |
| **Dependencies** | Requires `wide` + `bytemuck` | `fearless_simd` (zero external dependencies) |
| **Compiler Channel** | Stable / Nightly | Stable (MSRV 1.89) |

---

## Implementation Prototype Blueprint for `src/simd.rs`

```rust
use fearless_simd::{dispatch, Level, Simd};
use fearless_simd_macros::simd;

#[simd]
fn update_stats_simd_kernel<S: Simd>(
    simd: S,
    line: &[u8],
) -> (usize, usize, usize) {
    let mut gc_total = 0;
    let mut n_total = 0;
    let mut seq_chars_total = 0;

    let len = S::u8s::LEN;
    let mut chunks = line.chunks_exact(len);

    let v_case_mask = S::u8s::splat(simd, 0x20);
    let v_g = S::u8s::splat(simd, b'g');
    let v_c = S::u8s::splat(simd, b'c');
    let v_n = S::u8s::splat(simd, b'n');

    let v_space = S::u8s::splat(simd, b' ');
    let v_tab = S::u8s::splat(simd, b'\t');
    let v_nl = S::u8s::splat(simd, b'\n');
    let v_cr = S::u8s::splat(simd, b'\r');
    let v_dash = S::u8s::splat(simd, b'-');
    let v_dot = S::u8s::splat(simd, b'.');

    for chunk in chunks.by_ref() {
        let v_chunk = S::u8s::from_slice(simd, chunk);
        let v = v_chunk | v_case_mask;

        let is_gc = v.simd_eq(v_g) | v.simd_eq(v_c);
        let is_n = v.simd_eq(v_n);

        gc_total += is_gc.to_bitmask().count_ones() as usize;
        n_total += is_n.to_bitmask().count_ones() as usize;

        let is_skipped = v_chunk.simd_eq(v_space)
            | v_chunk.simd_eq(v_tab)
            | v_chunk.simd_eq(v_nl)
            | v_chunk.simd_eq(v_cr)
            | v_chunk.simd_eq(v_dash)
            | v_chunk.simd_eq(v_dot);

        let skipped_count = is_skipped.to_bitmask().count_ones() as usize;
        seq_chars_total += len - skipped_count;
    }

    let (gc, n, sc) = update_stats_scalar(chunks.remainder());
    (gc_total + gc, n_total + n, seq_chars_total + sc)
}

pub fn update_stats(line: &[u8], no_simd: bool) -> (usize, usize, usize) {
    if no_simd {
        return update_stats_scalar(line);
    }
    let level = Level::new();
    dispatch!(level, simd => update_stats_simd_kernel(simd, line))
}
```

---

## Conclusion & Recommendation

Migrating from `wide` to `fearless_simd` is a clean, drop-in refactoring of `src/simd.rs` that broadens CPU architecture optimization and removes `bytemuck` dependencies while preserving identical scalar verification semantics.
