# count-fasta-rs

This utility gives some stats on assembly reports

    Usage: count-fasta-rs [OPTIONS] [FASTA FILE (.fa .fasta .fna .zip .gz)]...

    Arguments:
      [FASTA FILE (.fa .fasta .fna .zip .gz)]...  

    Options:
      -c, --csv <CSV>              
      -d, --directory <DIRECTORY>  
      -h, --help                   Print help
      -V, --version                Print version

## Installation

Install the binary for your OS and architecture from the latest [release page](https://github.com/jovillarrealm/count-fasta-rs/releases) , or use the standalone installers. 

Or you can build from source with cargo.

    git clone https://github.com/jovillarrealm/count-fasta-rs
    cd count-fasta-rs
    cargo install --path .

## Usage 
When a csv file is not specified 

    count-fasta-rs ../cnsg-scripts/GENOMIC/GCA_024699835_Aphelenchoides-besseyi_AORJ.fna 

Output to stdout will look like this

```yaml
    Total length of sequence:       46759715 bp
    Total number of sequences:      32
    Average contig length is:       1461241 bp
    Largest contig:         18100598 bp
    Shortest contig:                214 bp
    N25 stats:                      25% of total sequence length is contained in the 1 sequences >= 18100598 bp
    N50 stats:                      50% of total sequence length is contained in the 2 sequences >= 16068654 bp
    N75 stats:                      75% of total sequence length is contained in the 3 sequences >= 10965501 bp
    Total GC count:                 19534458 bp
    GC %:                           41.78 %
    Number of Ns:                   2900
    Ns %:                           0.01 %
```

If the csv file was specified, then the created file will look like this
```rs
    filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage
    "GCA_024699835_Aphelenchoides-besseyi_AORJ.fna";46759715;32;1461241.09;18100598;214;16068654;41.78;2900;0.01
```

## Testing

To run the standard tests:

```bash
cargo test
```

### SIMD Cross-Architecture Testing

To verify the SIMD implementations (AVX2, AVX512, NEON) across different architectures using QEMU and `cross`:

1.  **Install `cross`**: `cargo install cross`
2.  **Ensure Docker/Podman** is running.
3.  **Run the cross-test script**:

```bash
./scripts/test_simd_cross.sh
```

This script uses `QEMU_CPU=max` to ensure the emulated CPUs support the required SIMD features. It tests both `aarch64-unknown-linux-gnu` (NEON) and `x86_64-unknown-linux-gnu` (AVX2/AVX512).

## Implementation deets

It uses [dist](https://github.com/axodotdev/cargo-dist) as a distribution tool for installation.

> [!IMPORTANT]
> To modify the installation path, it requires to modify install-path in dist-workspace.toml, and run ci again.
> `install-path = "CARGO_HOME"`

> [!NOTE]
> Known Bug:
> .gz files that were compressed with bgzip are read incorrectly.
