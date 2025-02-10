// kseq.h header (include at the top of your file)
#include "kseq.h"
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>

#define MAX_FILENAME 256
#define MAX_LINE 1000

typedef struct {
    char filename[MAX_FILENAME];
    size_t total_length;
    size_t sequence_count;
    size_t gc_count;
    size_t n_count;
    size_t n25;
    size_t n25_sequence_count;
    size_t n50;
    size_t n50_sequence_count;
    size_t n75;
    size_t n75_sequence_count;
    size_t largest_contig;
    size_t shortest_contig;
} AnalysisResults;

void init_results(AnalysisResults *results) {
    memset(results, 0, sizeof(AnalysisResults));
    results->shortest_contig = SIZE_MAX;
}

void process_sequence(const char *sequence, size_t length, AnalysisResults *results, size_t **lengths, size_t *lengths_count, size_t *lengths_capacity) {
    results->sequence_count++;
    results->total_length += length;
    
    for (size_t i = 0; i < length; i++) {
        char c = toupper(sequence[i]);
        if (c == 'G' || c == 'C') results->gc_count++;
        if (c == 'N') results->n_count++;
    }
    
    if (length > results->largest_contig) results->largest_contig = length;
    if (length < results->shortest_contig) results->shortest_contig = length;
    
    if (*lengths_count >= *lengths_capacity) {
        *lengths_capacity *= 2;
        *lengths = realloc(*lengths, *lengths_capacity * sizeof(size_t));
        if (!*lengths) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }
    (*lengths)[*lengths_count] = length;
    (*lengths_count)++;
}

int compare_size_t(const void *a, const void *b) {
    return (*(size_t*)b - *(size_t*)a);
}

void calc_nq_stats(size_t *lengths, size_t lengths_count, AnalysisResults *results) {
    qsort(lengths, lengths_count, sizeof(size_t), compare_size_t);
    
    size_t cumulative_length = 0;
    size_t cumulative_sequences = 0;
    
    for (size_t i = 0; i < lengths_count; i++) {
        cumulative_length += lengths[i];
        cumulative_sequences++;
        
        if (results->n25 == 0 && cumulative_length >= results->total_length / 4) {
            results->n25 = lengths[i];
            results->n25_sequence_count = cumulative_sequences;
        }
        if (results->n50 == 0 && cumulative_length >= results->total_length / 2) {
            results->n50 = lengths[i];
            results->n50_sequence_count = cumulative_sequences;
        }
        if (results->n75 == 0 && cumulative_length >= results->total_length * 3 / 4) {
            results->n75 = lengths[i];
            results->n75_sequence_count = cumulative_sequences;
            break;
        }
    }
}

void print_results(const AnalysisResults *results) {
    printf("\nTotal length of sequence:\t%zu bp\n", results->total_length);
    printf("Total number of sequences:\t%zu\n", results->sequence_count);
    printf("Average contig length is:\t%zu bp\n", results->total_length / results->sequence_count);
    printf("Largest contig:\t\t%zu bp\n", results->largest_contig);
    printf("Shortest contig:\t\t%zu bp\n", results->shortest_contig);
    printf("N25 stats:\t\t\t25%% of total sequence length is contained in the %zu sequences >= %zu bp\n",
           results->n25_sequence_count, results->n25);
    printf("N50 stats:\t\t\t50%% of total sequence length is contained in the %zu sequences >= %zu bp\n",
           results->n50_sequence_count, results->n50);
    printf("N75 stats:\t\t\t75%% of total sequence length is contained in the %zu sequences >= %zu bp\n",
           results->n75_sequence_count, results->n75);
    printf("Total GC count:\t\t\t%zu bp\n", results->gc_count);
    printf("GC %%:\t\t\t\t%.2f %%\n", (double)results->gc_count / results->total_length * 100.0);
    printf("Number of Ns:\t\t\t%zu\n", results->n_count);
    printf("Ns %%:\t\t\t\t%.2f %%\n", (double)results->n_count / results->total_length * 100.0);
}

void append_to_csv(const AnalysisResults *results, const char *csv_filename) {
    FILE *file = fopen(csv_filename, "a+");
    if (!file) {
        perror("Error opening CSV file");
        return;
    }

    // Check if the file is empty and write header if needed
    fseek(file, 0, SEEK_END);
    if (ftell(file) == 0) {
        fprintf(file, "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage\n");
    }

    // Lock the file
    struct flock lock;
    memset(&lock, 0, sizeof(lock));
    lock.l_type = F_WRLCK;
    fcntl(fileno(file), F_SETLKW, &lock);

    // Write data
    fprintf(file, "%s;%zu;%zu;%f;%zu;%zu;%zu;%f;%zu;%f\n",
            results->filename,
            results->total_length,
            results->sequence_count,
            (double)results->total_length / results->sequence_count,
            results->largest_contig,
            results->shortest_contig,
            results->n50,
            (double)results->gc_count / results->total_length * 100.0,
            results->n_count,
            (double)results->n_count / results->total_length * 100.0);

    // Unlock the file
    lock.l_type = F_UNLCK;
    fcntl(fileno(file), F_SETLK, &lock);

    fclose(file);
}

int main(int argc, char *argv[]) {
    char *csv_filename = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "c:")) != -1) {
        switch (opt) {
            case 'c':
                csv_filename = optarg;
                break;
            default:
                fprintf(stderr, "Usage: %s [-c csv_file] <fasta_file>\n", argv[0]);
                exit(1);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Expected FASTA file after options\n");
        exit(1);
    }

    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(argv[optind], "r");
    if (!fp) {
        perror("Error opening file");
        return 1;
    }
    seq = kseq_init(fp);

    AnalysisResults results;
    init_results(&results);
    strncpy(results.filename, basename(argv[optind]), MAX_FILENAME - 1);

    size_t *lengths = malloc(1000 * sizeof(size_t));
    size_t lengths_count = 0;
    size_t lengths_capacity = 1000;

    while ((l = kseq_read(seq)) >= 0) {
        process_sequence(seq->seq.s, seq->seq.l, &results, &lengths, &lengths_count, &lengths_capacity);
    }

    kseq_destroy(seq);
    gzclose(fp);

    calc_nq_stats(lengths, lengths_count, &results);
    
    if (csv_filename) {
        append_to_csv(&results, csv_filename);
    } else {
        print_results(&results);
    }

    free(lengths);
    return 0;
}

/*
Installation and Compilation Instructions:

1. Install zlib development files:
   On Ubuntu/Debian: sudo apt-get install zlib1g-dev
   On CentOS/RHEL: sudo yum install zlib-devel
   On macOS with Homebrew: brew install zlib

2. Download kseq.h:
   wget https://raw.githubusercontent.com/attractivechaos/klib/master/kseq.h

3. Compile the program:
   gcc -O2 -o count-fasta-c count-fasta-c.c -lz

4. Run the program:
   To print results to console: ./fasta_analyzer <path_to_fasta_file>
   To append results to CSV: ./fasta_analyzer -c output.csv <path_to_fasta_file>

Note: Make sure kseq.h is in the same directory as your C file when compiling.
*/
