#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "dna.h"


int main(int argc, char *argv[]) {
    if (argc == 2) {
        // Parse FASTA file
        parse_fasta(argv[1]);
    } else if (argc == 3) {
        // Decode SEQ to FASTA
        decode_seq_to_fasta(argv[1], argv[2]);
    } else {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s <input.fasta>           # Parse FASTA file\n", argv[0]);
        fprintf(stderr, "  %s <input.seq> <output.fasta>  # Convert SEQ to FASTA\n", argv[0]);
        return 1;
    }
    return 0;
}