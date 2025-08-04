#ifndef TWO_BIT_DNA
#define TWO_BIT_DNA

#include <stdint.h>

// Constants
#define SEQ_SIGNATURE "SEQ\x01"
static const uint8_t VERSION = 0x01;

// Types
#define TYPE_DNA 1
#define TYPE_RNA 2

// Function Prototypes

uint8_t base_to_bits(char base);
int detect_rna(const char *seq);
void write_seq_file(const char *id, const char *desc, const char *seq);
void parse_fasta(const char *filename);

char bits_to_base(uint8_t bits, uint8_t type);
void decode_seq_to_fasta(const char *seqfile, const char *fastafile);


#endif
