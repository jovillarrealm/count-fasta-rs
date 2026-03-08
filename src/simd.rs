use std::simd::prelude::*;

const LOOKUP: [u8; 256] = {
    let mut table = [0u8; 256];
    // Bit 0: GC
    // Bit 1: N
    // Bit 2: Skip (whitespace, gaps)
    
    // GC
    table[b'G' as usize] = 1;
    table[b'g' as usize] = 1;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    
    // N
    table[b'N' as usize] = 2;
    table[b'n' as usize] = 2;
    
    // Skip (whitespace and gaps)
    table[b' ' as usize] = 4;
    table[b'\t' as usize] = 4;
    table[b'\n' as usize] = 4;
    table[b'\r' as usize] = 4;
    table[b'-' as usize] = 4;
    table[b'.' as usize] = 4;
    
    table
};

pub fn update_stats(line: &[u8], no_simd: bool) -> (usize, usize, usize) {
    if no_simd {
        return update_stats_scalar(line);
    }

    update_stats_simd(line)
}

fn update_stats_scalar(line: &[u8]) -> (usize, usize, usize) {
    let mut gc = 0;
    let mut n = 0;
    let mut seq_chars = 0;
    for &b in line {
        let val = LOOKUP[b as usize];
        gc += (val & 1) as usize;
        n += ((val & 2) >> 1) as usize;
        seq_chars += (1 - (val >> 2)) as usize;
    }
    (gc, n, seq_chars)
}

fn update_stats_simd(line: &[u8]) -> (usize, usize, usize) {
    let mut gc_total = 0;
    let mut n_total = 0;
    let mut seq_chars_total = 0;

    let (head, mid, tail) = line.as_simd::<32>();

    // Process head
    let (gc, n, sc) = update_stats_scalar(head);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    let v_case_mask = u8x32::splat(0x20);
    let v_g = u8x32::splat(b'g');
    let v_c = u8x32::splat(b'c');
    let v_n = u8x32::splat(b'n');
    
    let v_space = u8x32::splat(b' ');
    let v_tab = u8x32::splat(b'\t');
    let v_nl = u8x32::splat(b'\n');
    let v_cr = u8x32::splat(b'\r');
    let v_dash = u8x32::splat(b'-');
    let v_dot = u8x32::splat(b'.');

    for &chunk in mid {
        let v = chunk | v_case_mask;
        
        let is_g = v.simd_eq(v_g);
        let is_c = v.simd_eq(v_c);
        let is_n = v.simd_eq(v_n);
        
        let is_gc = is_g | is_c;
        gc_total += is_gc.to_bitmask().count_ones() as usize;
        n_total += is_n.to_bitmask().count_ones() as usize;

        // Count skipped
        let s1 = chunk.simd_eq(v_space);
        let s2 = chunk.simd_eq(v_tab);
        let s3 = chunk.simd_eq(v_nl);
        let s4 = chunk.simd_eq(v_cr);
        let s5 = chunk.simd_eq(v_dash);
        let s6 = chunk.simd_eq(v_dot);
        
        let is_skipped = s1 | s2 | s3 | s4 | s5 | s6;
        let skipped_count = is_skipped.to_bitmask().count_ones() as usize;
        seq_chars_total += 32 - skipped_count;
    }

    // Process tail
    let (gc, n, sc) = update_stats_scalar(tail);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    (gc_total, n_total, seq_chars_total)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Simple Xorshift RNG for testing without extra dependencies
    struct SimpleRng {
        state: u64,
    }

    impl SimpleRng {
        fn new(seed: u64) -> Self {
            Self { state: seed }
        }

        fn next_u8(&mut self) -> u8 {
            self.state ^= self.state << 13;
            self.state ^= self.state >> 7;
            self.state ^= self.state << 17;
            (self.state & 0xFF) as u8
        }

        fn fill_bytes(&mut self, buf: &mut [u8]) {
            for b in buf.iter_mut() {
                *b = self.next_u8();
            }
        }
    }

    fn check_consistency(input: &[u8]) {
        let scalar_res = update_stats_scalar(input);
        // Test normal dispatch
        let dispatch_res = update_stats(input, false);
        assert_eq!(scalar_res, dispatch_res, "Dispatched (SIMD enabled) result should match scalar result for len {}", input.len());
        
        // Test no_simd flag
        let no_simd_res = update_stats(input, true);
        assert_eq!(scalar_res, no_simd_res, "Dispatched (SIMD disabled) result should match scalar result for len {}", input.len());
    }

    #[test]
    fn test_update_stats_consistency_basic() {
        let patterns: Vec<&[u8]> = vec![
            b"AaCcGgTtNnRrYyWwSsMmKkHhBbVvDd-. \t\n\r",
            b"G",
            b"GC",
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // 32 bytes
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // 33 bytes
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",  // 31 bytes
            b"",
        ];

        for input in patterns {
            check_consistency(input);
        }
    }

    #[test]
    fn test_fuzz_update_stats() {
        let mut rng = SimpleRng::new(12345);
        for _ in 0..100 {
            // Generate random length up to ~1MB
            // (0..255 * 4096) + (0..255 * 16) + (0..255)
            let len = (rng.next_u8() as usize) * 4096 
                    + (rng.next_u8() as usize) * 16 
                    + (rng.next_u8() as usize); 
            
            let mut buf = vec![0u8; len];
            rng.fill_bytes(&mut buf);
            
            // Constrain bytes to likely FASTA characters to test specific logic paths
            for b in buf.iter_mut() {
                let r = *b % 20;
                *b = match r {
                    0..=3 => b'G',
                    4..=7 => b'C', 
                    8..=9 => b'N',
                    10 => b'n',
                    11 => b'g',
                    12 => b'c',
                    13 => b'\n',
                    14 => b' ',
                    15 => b'-',
                    _ => b'A', // Other chars
                };
            }

            check_consistency(&buf);
        }
    }

    #[test]
    fn test_alignment_and_offsets() {
        // Allocate a larger buffer to play with alignments and larger chunks
        let mut base_buf = vec![b'A'; 4096 * 4]; 
        // Fill with some pattern
        for (i, b) in base_buf.iter_mut().enumerate() {
            if i % 3 == 0 { *b = b'G'; }
            if i % 5 == 0 { *b = b'N'; }
            if i % 7 == 0 { *b = b'\n'; }
        }

        // Test various sub-slices to force different alignments and head/tail combinations
        for start in 0..64 {
            for len in 0..512 {
                if start + len > base_buf.len() { break; }
                let slice = &base_buf[start..start+len];
                check_consistency(slice);
            }
        }
    }

    #[test]
    fn test_vector_boundaries() {
        let max_vec_size = 64;
        let base_buf = vec![b'G'; max_vec_size * 3];

        // Test lengths around vector sizes
        for size in [16, 32, 64] {
            if size * 2 > base_buf.len() { continue; }
            
            check_consistency(&base_buf[0..size]);        // Exact size
            check_consistency(&base_buf[0..size - 1]);    // Size - 1
            check_consistency(&base_buf[0..size + 1]);    // Size + 1
            check_consistency(&base_buf[0..size * 2]);    // Double size
            check_consistency(&base_buf[0..size * 2 - 1]);
            check_consistency(&base_buf[0..size * 2 + 1]);
        }
    }

    #[test]
    fn test_resilience_random_junk() {
        let mut rng = SimpleRng::new(9999);
        for _ in 0..100 {
            let len = (rng.next_u8() as usize) * 64 + (rng.next_u8() as usize);
            let mut buf = vec![0u8; len];
            rng.fill_bytes(&mut buf);
            check_consistency(&buf);
        }
    }

    #[test]
    fn test_resilience_edge_cases() {
        let edge_cases: Vec<&[u8]> = vec![
            b"ACGT\0ACGT", 
            b"\0\0\0",
            b"ACGT\xFF\x80ACGT",
            b"\x01\x02\x03\n\r\t\x0B\x0C",
            b">seq1\n>seq2\n>seq3",
            b"   \t\t\n\n\r\r",
            b"AgC\xFFT\x00N-",
        ];

        for case in edge_cases {
            check_consistency(case);
        }
    }
}
