# xumi

Extract subsequences from read alignments at specified genomic regions,
with full CIGAR awareness. Originally developed for inline UMI extraction.

## Overview

**xumi** reads a SAM/BAM/CRAM file and, for each mapped read, extracts the
subsequence(s) corresponding to one or more user-specified reference regions.
It handles complex CIGAR strings correctly — insertions, deletions, soft clips,
and all standard alignment operations are taken into account.

## Motivation

Existing tools solve adjacent problems but not this specific one:

| Tool | What it does | Gap |
|---|---|---|
| `samtools view` | Extracts whole reads overlapping a region | Does not isolate the subsequence within the region |
| `bedtools getfasta` | Extracts the reference sequence at a region | Returns the reference, not per-read sequences |
| UMI-tools | Extracts UMIs from fixed read positions (first/last N bases) | Cannot handle UMIs at arbitrary genomic coordinates |
| `samtools consensus` | Builds a consensus over a region | Collapses reads; no per-read output |

**xumi** fills this gap: given a mapped read and a reference region, return
the portion of *that read* corresponding to *that region*, with full
CIGAR awareness.

### Key features

- Full CIGAR awareness (M, =, X, I, D, N, S, H, P)
- Multiple input regions per run
- Control over boundary insertion inclusion
- Aligned-bases-only mode (no insertions)
- TSV and FASTA output, in wide or long layout
- Reads from stdin, writes to stdout — pipeable with samtools

## Installation

### With pip (in a virtual environment)

```bash
python -m venv .venv
source .venv/bin/activate
pip install .
```

### With conda

```bash
conda create -n xumi python=3.10
conda activate xumi
pip install .
```

### Requirements

- Python ≥ 3.10 (the code itself does not use features newer than 3.7,
  but it has only been tested on 3.10)
- [pysam](https://github.com/pysam-developers/pysam) ≥ 0.20

## Quick start

```bash
# Single region, BAM input
xumi -r chr1:100-200 mapped.bam

# Multiple regions
xumi -r chr1:100-200,chr1:1000-1100 mapped.bam

# Regions from a BED file
xumi -R regions.bed mapped.bam

# Pipe from samtools
samtools view -u mapped.bam chr1 | xumi -r chr1:100-200
```

## Usage

```text
xumi [-h] [-V] [-H] [-r REGIONS] [-R in.bed]
     [-a] [-b {both,left,right,none}]
     [-O {tsv,tsv-long,fasta,fasta-long}]
     [-o FILE] [aln.sam|aln.bam|aln.cram]
```

### Options

| Option | Description |
|---|---|
| `-r`, `--regions` | Comma-separated regions as `RNAME:START-END` (1-based) |
| `-R`, `--regions-file` | BED file containing regions |
| `-a`, `--aligned-only` | Only output reference-aligned bases (M/=/X) |
| `-b`, `--boundary-insertions` | Include boundary insertions: `both` (default), `left`, `right`, `none` |
| `-O`, `--output-format` | Output format (see below) |
| `-o`, `--output` | Output file (default: stdout, `.gz` for gzipped) |
| `-H`, `--no-header` | Suppress header line(s) |
| `-V`, `--version` | Show version |

At least one of `-r` or `-R` is required.

### Output formats

| Format | Layout | Description |
|---|---|---|
| `tsv` (default) | wide | One row per read, one column per region |
| `tsv-long` | long | One row per read × region pair |
| `fasta` | wide | One record per read, regions joined with `-` |
| `fasta-long` | long | One record per read × region pair |

Reads with no extracted sequence are omitted.

## Extraction modes

### Default: query slice

Returns the contiguous portion of the read spanning the region, including
any insertions between the first and last aligned bases. Boundary
insertions (insertions immediately adjacent to the region edges) are
controlled with `--boundary-insertions`.

### Aligned-only (`--aligned-only`)

Returns only the reference-aligned bases (M/=/X operations) overlapping
the region. Insertions are excluded entirely. Useful when you need the
exact reference-matching content.

### Boundary insertions (`--boundary-insertions`)

When extracting contiguous regions, insertions at the boundaries can
belong to either side. Use this option to control which region gets them:

```bash
# Extract two adjacent regions without double-counting boundary insertions
xumi -r chr1:10-19,chr1:20-29 regions.bed --boundary-insertions left mapped.bam
```

| Value | Left boundary | Right boundary |
|---|---|---|
| `both` (default) | included | included |
| `left` | included | excluded |
| `right` | excluded | included |
| `none` | excluded | excluded |

## Examples

### Extract UMIs from known positions

```bash
# UMI at chr1:1-8
xumi -r chr1:1-8 -O tsv mapped.bam > umis.tsv
```

### Multiple regions to TSV

```bash
xumi -r chr1:1-8,chr1:101-108 mapped.bam
```

Output:
```text
#qname	chr1:1-8	chr1:101-108
read1	ACGTACGT	TGCATGCA
read2	ACGTACG	TGCATGCAT
```

### FASTA output for downstream tools

```bash
xumi -r chr1:1-8 -O fasta-long mapped.bam > extracted.fa
```

Output:
```
>read1::chr1:1-8
ACGTACGT
>read2::chr1:1-8
ACGTACG
```

## Supplementary and secondary alignments

Supplementary and secondary alignments are processed like any other
mapped read. A read name may therefore appear multiple times in the
output when a read has multiple alignments.

## Debugging

Set the `XUMI_DEBUG` environment variable for full tracebacks on error:

```bash
XUMI_DEBUG=1 xumi -r chr1:100-200 mapped.bam
```

## License

[MIT](LICENSE)

## Author

Rafael NAVAZA — Institut Pasteur, Plateforme de Cristallographie

rnavaza@pasteur.fr
