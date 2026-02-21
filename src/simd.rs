
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

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

pub fn update_stats(line: &[u8]) -> (usize, usize, usize) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx512bw") {
            return unsafe { update_stats_avx512(line) };
        }
        if is_x86_feature_detected!("avx2") {
            return unsafe { update_stats_avx2(line) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            return unsafe { update_stats_neon(line) };
        }
    }
    
    update_stats_scalar(line)
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

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn update_stats_avx2(line: &[u8]) -> (usize, usize, usize) {
    let mut gc_total = 0;
    let mut n_total = 0;
    let mut seq_chars_total = 0;

    let (head, mid, tail) = unsafe { line.align_to::<__m256i>() };

    // Process head
    let (gc, n, sc) = update_stats_scalar(head);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    let v_case_mask = _mm256_set1_epi8(0x20);
    let v_g = _mm256_set1_epi8(b'g' as i8);
    let v_c = _mm256_set1_epi8(b'c' as i8);
    let v_n = _mm256_set1_epi8(b'n' as i8);
    
    let v_space = _mm256_set1_epi8(b' ' as i8);
    let v_tab = _mm256_set1_epi8(b'\t' as i8);
    let v_nl = _mm256_set1_epi8(b'\n' as i8);
    let v_cr = _mm256_set1_epi8(b'\r' as i8);
    let v_dash = _mm256_set1_epi8(b'-' as i8);
    let v_dot = _mm256_set1_epi8(b'.' as i8);

    for &chunk in mid {
        let v = _mm256_or_si256(chunk, v_case_mask);
        
        let is_g = _mm256_cmpeq_epi8(v, v_g);
        let is_c = _mm256_cmpeq_epi8(v, v_c);
        let is_n = _mm256_cmpeq_epi8(v, v_n);
        
        let is_gc = _mm256_or_si256(is_g, is_c);
        gc_total += _mm256_movemask_epi8(is_gc).count_ones() as usize;
        n_total += _mm256_movemask_epi8(is_n).count_ones() as usize;

        // Count skipped
        let s1 = _mm256_cmpeq_epi8(chunk, v_space);
        let s2 = _mm256_cmpeq_epi8(chunk, v_tab);
        let s3 = _mm256_cmpeq_epi8(chunk, v_nl);
        let s4 = _mm256_cmpeq_epi8(chunk, v_cr);
        let s5 = _mm256_cmpeq_epi8(chunk, v_dash);
        let s6 = _mm256_cmpeq_epi8(chunk, v_dot);
        
        let is_skipped = _mm256_or_si256(
            _mm256_or_si256(_mm256_or_si256(s1, s2), _mm256_or_si256(s3, s4)),
            _mm256_or_si256(s5, s6)
        );
        let skipped_count = _mm256_movemask_epi8(is_skipped).count_ones() as usize;
        seq_chars_total += 32 - skipped_count;
    }

    // Process tail
    let (gc, n, sc) = update_stats_scalar(tail);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    (gc_total, n_total, seq_chars_total)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx512bw")]
unsafe fn update_stats_avx512(line: &[u8]) -> (usize, usize, usize) {
    let mut gc_total = 0;
    let mut n_total = 0;
    let mut seq_chars_total = 0;

    let (head, mid, tail) = unsafe { line.align_to::<__m512i>() };

    // Process head
    let (gc, n, sc) = update_stats_scalar(head);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    let v_case_mask = _mm512_set1_epi8(0x20);
    let v_g = _mm512_set1_epi8(b'g' as i8);
    let v_c = _mm512_set1_epi8(b'c' as i8);
    let v_n = _mm512_set1_epi8(b'n' as i8);
    
    let v_space = _mm512_set1_epi8(b' ' as i8);
    let v_tab = _mm512_set1_epi8(b'\t' as i8);
    let v_nl = _mm512_set1_epi8(b'\n' as i8);
    let v_cr = _mm512_set1_epi8(b'\r' as i8);
    let v_dash = _mm512_set1_epi8(b'-' as i8);
    let v_dot = _mm512_set1_epi8(b'.' as i8);

    for &chunk in mid {
        let v = _mm512_or_si512(chunk, v_case_mask);
        
        let is_g = _mm512_cmpeq_epi8_mask(v, v_g);
        let is_c = _mm512_cmpeq_epi8_mask(v, v_c);
        let is_n = _mm512_cmpeq_epi8_mask(v, v_n);
        
        let is_gc = is_g | is_c;
        gc_total += is_gc.count_ones() as usize;
        n_total += is_n.count_ones() as usize;

        // Count skipped
        let s1 = _mm512_cmpeq_epi8_mask(chunk, v_space);
        let s2 = _mm512_cmpeq_epi8_mask(chunk, v_tab);
        let s3 = _mm512_cmpeq_epi8_mask(chunk, v_nl);
        let s4 = _mm512_cmpeq_epi8_mask(chunk, v_cr);
        let s5 = _mm512_cmpeq_epi8_mask(chunk, v_dash);
        let s6 = _mm512_cmpeq_epi8_mask(chunk, v_dot);
        
        let is_skipped = s1 | s2 | s3 | s4 | s5 | s6;
        let skipped_count = is_skipped.count_ones() as usize;
        seq_chars_total += 64 - skipped_count;
    }

    // Process tail
    let (gc, n, sc) = update_stats_scalar(tail);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    (gc_total, n_total, seq_chars_total)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn update_stats_neon(line: &[u8]) -> (usize, usize, usize) {
    let mut gc_total = 0;
    let mut n_total = 0;
    let mut seq_chars_total = 0;

    let (head, mid, tail) = unsafe { line.align_to::<uint8x16_t>() };

    // Process head
    let (gc, n, sc) = update_stats_scalar(head);
    gc_total += gc;
    n_total += n;
    seq_chars_total += sc;

    let v_case_mask = vdupq_n_u8(0x20);
    let v_g = vdupq_n_u8(b'g');
    let v_c = vdupq_n_u8(b'c');
    let v_n = vdupq_n_u8(b'n');
    
    let v_space = vdupq_n_u8(b' ');
    let v_tab = vdupq_n_u8(b'\t');
    let v_nl = vdupq_n_u8(b'\n');
    let v_cr = vdupq_n_u8(b'\r');
    let v_dash = vdupq_n_u8(b'-');
    let v_dot = vdupq_n_u8(b'.');

    for &chunk in mid {
        let v = vorrq_u8(chunk, v_case_mask);
        
        let is_g = vceqq_u8(v, v_g);
        let is_c = vceqq_u8(v, v_c);
        let is_n = vceqq_u8(v, v_n);
        
        let is_gc = vorrq_u8(is_g, is_c);
        
        let v_one = vdupq_n_u8(1);
        let gc_counts = vandq_u8(is_gc, v_one);
        let n_counts = vandq_u8(is_n, v_one);
        
        gc_total += vaddlvq_u8(gc_counts) as usize;
        n_total += vaddlvq_u8(n_counts) as usize;

        // Count skipped
        let s1 = vceqq_u8(chunk, v_space);
        let s2 = vceqq_u8(chunk, v_tab);
        let s3 = vceqq_u8(chunk, v_nl);
        let s4 = vceqq_u8(chunk, v_cr);
        let s5 = vceqq_u8(chunk, v_dash);
        let s6 = vceqq_u8(chunk, v_dot);
        
        let is_skipped = vorrq_u8(
            vorrq_u8(vorrq_u8(s1, s2), vorrq_u8(s3, s4)),
            vorrq_u8(s5, s6)
        );
        let skipped_counts = vandq_u8(is_skipped, v_one);
        let skipped_count = vaddlvq_u8(skipped_counts) as usize;
        
        seq_chars_total += 16 - skipped_count;
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
        let dispatch_res = update_stats(input);
        
        assert_eq!(scalar_res, dispatch_res, "Dispatched result should match scalar result for len {}", input.len());
        
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if is_x86_feature_detected!("avx2") {
            let avx2_res = unsafe { update_stats_avx2(input) };
            assert_eq!(scalar_res, avx2_res, "AVX2 result should match scalar result for len {}", input.len());
        }

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if is_x86_feature_detected!("avx512bw") {
             let avx512_res = unsafe { update_stats_avx512(input) };
             assert_eq!(scalar_res, avx512_res, "AVX512 result should match scalar result for len {}", input.len());
        }

        #[cfg(target_arch = "aarch64")]
        if std::arch::is_aarch64_feature_detected!("neon") {
            let neon_res = unsafe { update_stats_neon(input) };
            assert_eq!(scalar_res, neon_res, "NEON result should match scalar result for len {}", input.len());
        }
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
        // align_to depends on the pointer address, so iterating offsets is crucial.
        for start in 0..64 {
            // Check lengths up to 512 bytes to cover multiple vector widths
            for len in 0..512 {
                if start + len > base_buf.len() { break; }
                let slice = &base_buf[start..start+len];
                check_consistency(slice);
            }
        }
    }

    #[test]
    fn test_vector_boundaries() {
        let max_vec_size = 64; // AVX512 is 64 bytes
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
        // Test with completely random junk (0-255) to ensure no crashes on invalid ASCII/UTF-8
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
            // Embedded nulls
            b"ACGT\0ACGT", 
            b"\0\0\0",
            // High bit set (extended ASCII / invalid UTF-8 start bytes)
            b"ACGT\xFF\x80ACGT",
            // All control characters
            b"\x01\x02\x03\n\r\t\x0B\x0C",
            // "Oops all headers" like pattern (though logic doesn't parse headers here, just stats)
            b">seq1\n>seq2\n>seq3",
            // Just whitespace
            b"   \t\t\n\n\r\r",
            // Alternating valid/invalid tightly
            b"AgC\xFFT\x00N-",
        ];

        for case in edge_cases {
            check_consistency(case);
        }
    }
}
