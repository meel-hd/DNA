# Optimized DNA/RNA Storage

## Sequence Encoder and Decoder

This project provides tools to:

* Convert multi-sequence FASTA files into individual compressed `.seq` binary files.
* Encode DNA/RNA sequences using 2-bit representations.
* Preserve metadata, type (DNA/RNA), and support versioning.
* Reconstruct `.fasta` files from `.seq` binaries.

## Format Specification

### `.seq` File Structure

```
Offset  Size      Description
0       4         Signature: "SEQ\x01"
4       1         Version (currently 0x01)
5       2         Metadata length (uint32_t)
9       N         Metadata string (FASTA header)
9+N     1         Type: 1 = DNA, 2 = RNA
10+N    4         Sequence length in bases (uint64_t)
18+N    ?         Sequence data (2 bits per base, packed)
```

Bases are encoded as follows:

| Base | Bits | Decimal |
| ---- | ---- | ------- |
| A    | 00   |    0    |
| C    | 01   |    1    |
| G    | 10   |    2    |
| T/U  | 11   |    3    |

### Encoding Notes

* RNA sequences are detected if they contain `U` and converted to type 2.
* The encoder replaces `U` with `T` internally for 2-bit packing.
* Metadata is the entire FASTA header line (starting with `>`).

## Usage


This program supports two modes based on the number of arguments:

```sh
# Converts FASTA to .seq files
./program <input.fasta> 
# Converts a .seq file back to FASTA
./program <input.seq> <output.fasta>  
```
- If only one argument is given, it treats it as a FASTA file and parses it into .seq files.
- If two arguments are given, it decodes the .seq file into a valid FASTA file.

## Example FASTA File

```
>example_sequence_53
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
```

## Compression Efficiency at Scale

The following demonstrates projected compression results when scaling to a large dataset:

- **Original FASTA file size:** 39.9 GB  
- **After `.seq` encoding:** 9.8 GB  
- **After 7zip compression:** 8.2 GB

This shows that with 2-bit encoding and efficient archiving:
- The encoded file is ~75% smaller than the original FASTA.
- The archived `.dna` file further reduces storage requirements.
- Metadata remains accessible without decompressing the entire archive.

This scale makes the format practical for handling large genomic datasets such as entire chromosomes or transcriptome libraries.


## Future Features

* Combine multiple `.seq` files into a single `.dna` archive/
* Add `.dna` archive reader/unpacker.
* Optional compression for `.dna` archives with metadata.
* GUI Tools.

**Note:**  `7zip` is a good fit for its robust archiving capabilities, especially its support for bundling both binary `.seq` files and accompanying metadata (e.g., JSON) to a `.dna` archive. This allows for efficient storage and distribution, while also enabling users,gui, or tools to quickly extract metadata files without unpacking the entire (potentially large) archive.

---

## Handling Biological Data Quirks

FASTA files from real-world sources often contain edge cases that must be handled with care. Our format currently uses a compact 2-bit encoding scheme for DNA/RNA sequences, which assumes only four standard bases: `A`, `C`, `G`, and `T`/`U`.

However, actual sequence data frequently includes the following quirks:

### 1. `N` Bases (Unknown/Any)
- `N` is a standard IUPAC nucleotide representing any base (A/C/G/T).
- Common in genome gaps, unresolved regions, and low-quality reads.
- **Solution:** We store the positions of `N` bases in the metadata to preserve full sequence information.

### 2. Other IUPAC Ambiguity Codes
- Codes like `R`, `Y`, `S`, `W`, etc., represent multiple possible bases.
- These are less common but appear in some variant-rich datasets.
- **Current Handling:** These bases are **not yet supported** and will cause the encoder to skip or abort. A warning is logged.
- **Future Plan:** Consider expanding support or filtering them pre-encoding.

### 3. Mixed or Lowercase Bases: Done ✅
- Some FASTA files use lowercase letters to represent soft-masked regions (e.g., repeats).
- **Current Handling:** All bases are normalized to uppercase.

### 4. Sequence Wrapping and Line Endings:  Done ✅
- Lines may be wrapped at different lengths (60, 80, or unwrapped).
- Files may include Windows (`\r\n`) or UNIX (`\n`) line endings.
- **Current Handling:** The parser reads sequences continuously, ignoring line breaks.

### 5. Multiple Sequences per FASTA
- Many FASTA files contain multiple records.
- **Handling:** Each sequence is extracted and encoded into a separate `.seq` file with metadata and type information.

---

**Note:** Our `.seq` encoder is intentionally strict to ensure reliable 2-bit packing. However, these quirks must be resolved (filtered, tracked, or logged) to avoid losing information or causing encoding failures.


## Build

Compile with a C compiler:

```sh
make
```

---

## License

Apache 2.0 [License](LICENSE)
