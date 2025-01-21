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
If you are setting a -d directory, you should probably set a -c csv file to go with it

## Implementation deets

It uses tokio to handle I/O asynchronously. It uses flate2 to handle gzipped files.
