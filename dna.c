// dna.c --- DNA Sequence Processing Functions
// This file contains functions for converting DNA sequences
// between the custom binary format and FASTA format.
// It includes functions for reading, writing, and parsing sequences,
// as well as converting between base representations.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "dna.h"


/**
 * Converts a single base character to its 2-bit representation.
 * @param base The base character ('A', 'C', 'G', 'T', 'U').
 * @return The 2-bit representation of the base.
 * @note 'U' is treated as 'T' for DNA sequences.
 *       If an invalid base is provided, the program will exit with an error.
 */
uint8_t base_to_bits(char base) {
    switch (base) {
        case 'A': return 0b00;
        case 'C': return 0b01;
        case 'G': return 0b10;
        case 'T': case 'U': return 0b11;
        default:
            fprintf(stderr, "Invalid base: %c\n", base);
            exit(EXIT_FAILURE);
    }
}

/**
 * Detects if the sequence is RNA based on the presence of 'U'.
 * @param seq The sequence string.
 * @return 1 if RNA, 0 if DNA.
 */
int detect_rna(const char *seq) {
    for (int i = 0; seq[i]; i++) {
        if (toupper(seq[i]) == 'U') return 1;
    }
    return 0;
}

/**
 * Writes a sequence to a .seq file in the custom binary format.
 * @param id The identifier for the sequence.
 * @param desc The description of the sequence.
 * @param seq The sequence string.
 */
void write_seq_file(const char *id, const char *desc, const char *seq) {
    char filename[256];
    snprintf(filename, sizeof(filename), "%s.seq", id);

    FILE *out = fopen(filename, "wb");
    if (!out) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }

    // Normalize sequence
    size_t len = strlen(seq);
    char *cleaned = malloc(len + 1);
    int j = 0;
    for (int i = 0; seq[i]; i++) {
        char c = toupper(seq[i]);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U') {
            cleaned[j++] = (c == 'U') ? 'T' : c;
        }
    }
    cleaned[j] = '\0';

    uint64_t seq_len = j;
    uint32_t meta_len = (uint32_t)strlen(desc);
    uint8_t type = detect_rna(seq) ? TYPE_RNA : TYPE_DNA;

    // Write signature and version
    fwrite(SEQ_SIGNATURE, 1, 4, out);
    fwrite(&VERSION, 1, 1, out);

    // Metadata
    fwrite(&meta_len, sizeof(uint32_t), 1, out);
    fwrite(desc, 1, meta_len, out);

    // Type and length
    fwrite(&type, 1, 1, out);
    fwrite(&seq_len, sizeof(uint64_t), 1, out);

    // Sequence data
    uint8_t buffer = 0;
    int bits = 0;
    for (uint64_t i = 0; i < seq_len; i++) {
        buffer <<= 2;
        buffer |= base_to_bits(cleaned[i]);
        bits += 2;
        if (bits == 8) {
            fwrite(&buffer, 1, 1, out);
            buffer = 0;
            bits = 0;
        }
    }
    if (bits > 0) {
        buffer <<= (8 - bits);
        fwrite(&buffer, 1, 1, out);
    }

    fclose(out);
    free(cleaned);
}

/**
 * Parses a FASTA file and writes each sequence to a .seq file.
 * @param filename The name of the FASTA file.
 */
void parse_fasta(const char *filename) {
    FILE *in = fopen(filename, "r");
    if (!in) {
        perror("Failed to open FASTA file");
        exit(EXIT_FAILURE);
    }

    char line[4096], *seq = NULL, *desc = NULL, *id = NULL;
    size_t seqcap = 0, seqlen = 0;

    while (fgets(line, sizeof(line), in)) {
        if (line[0] == '>') {
            if (id && desc && seq) {
                seq[seqlen] = '\0';
                write_seq_file(id, desc, seq);
                free(seq); seq = NULL;
                free(id); id = NULL;
                free(desc); desc = NULL;
                seqcap = seqlen = 0;
            }

            line[strcspn(line, "\r\n")] = 0;
            desc = strdup(line);
            char *space = strchr(desc, ' ');
            if (space) *space = '\0';
            id = strdup(desc + 1); // skip '>'
            if (space) *space = ' ';
        } else {
            size_t len = strlen(line);
            while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) len--;
            if (seqlen + len + 1 > seqcap) {
                seqcap = (seqlen + len + 1) * 2;
                seq = realloc(seq, seqcap);
            }
            for (size_t i = 0; i < len; i++) seq[seqlen++] = line[i];
        }
    }

    if (id && desc && seq) {
        seq[seqlen] = '\0';
        write_seq_file(id, desc, seq);
        free(seq); free(id); free(desc);
    }

    fclose(in);
}


/**
 * Converts a 2-bit representation back to a base character.
 * @param bits The 2-bit representation of the base.
 * @param type The type of sequence (TYPE_DNA or TYPE_RNA).
 * @return The base character ('A', 'C', 'G', 'T', 'U').
 * @note 'U' is returned for RNA sequences, 'T' for DNA sequences.
 */
char bits_to_base(uint8_t bits, uint8_t type) {
    switch (bits & 0b11) {
        case 0b00: return 'A';
        case 0b01: return 'C';
        case 0b10: return 'G';
        case 0b11: return (type == TYPE_RNA) ? 'U' : 'T';
        default:
            fprintf(stderr, "Invalid bits: %u\n", bits);
            exit(EXIT_FAILURE);
    }
}

/**
 * Decodes a .seq file back into a FASTA file.
 * @param seqfile The input .seq file.
 * @param fastafile The output FASTA file.
 */
void decode_seq_to_fasta(const char *seqfile, const char *fastafile) {
    FILE *in = fopen(seqfile, "rb");
    if (!in) {
        perror("Failed to open .seq file");
        exit(EXIT_FAILURE);
    }

    char signature[5] = {0};
    fread(signature, 1, 4, in);
    if (memcmp(signature, SEQ_SIGNATURE, 4) != 0) {
        fprintf(stderr, "Invalid .seq signature\n");
        exit(EXIT_FAILURE);
    }

    uint8_t version;
    fread(&version, 1, 1, in);

    uint32_t meta_len;
    fread(&meta_len, sizeof(uint32_t), 1, in);

    char *meta = malloc(meta_len + 1);
    fread(meta, 1, meta_len, in);
    meta[meta_len] = '\0';

    uint8_t type;
    fread(&type, 1, 1, in);

    uint64_t seqlen;
    fread(&seqlen, sizeof(uint64_t), 1, in);

    FILE *out = fopen(fastafile, "w");
    if (!out) {
        perror("Failed to open output FASTA file");
        exit(EXIT_FAILURE);
    }

    fprintf(out, "%s\n", meta);

    int bases_written = 0;
    uint8_t byte;
    while (fread(&byte, 1, 1, in) == 1 && bases_written < seqlen) {
        for (int shift = 6; shift >= 0 && bases_written < seqlen; shift -= 2) {
            char base = bits_to_base((byte >> shift) & 0b11, type);
            fputc(base, out);
            bases_written++;
            if (bases_written % 60 == 0) fputc('\n', out);
        }
    }
    if (bases_written % 60 != 0) fputc('\n', out);

    fclose(in);
    fclose(out);
    free(meta);
}
