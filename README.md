# xumi

![Python](https://img.shields.io/badge/python-3.10-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![GitHub](https://img.shields.io/github/release/fravadona/xumi)
[![DOI](https://zenodo.org/badge/1172913507.svg)](https://doi.org/10.5281/zenodo.18905757)

Extract subsequences from read alignments at specified genomic regions,
with full CIGAR awareness.

## Overview

Originally developed for inline dual-UMI extraction, but applicable to any
region-based sequence extraction from aligned reads.  

**xumi** reads a SAM/BAM/CRAM file and, for each mapped read, extracts the
subsequence(s) corresponding to one or more user-specified reference regions
using the CIGAR information contained in the read record.

### Motivation

**xumi** fills a gap in the ecosystem, as existing tools solve adjacent problems:

| Tool | Limitation |
|---|---|
| `samtools view` | Returns whole reads, not the subsequence within a region |
| `bedtools getfasta` | Returns the reference sequence, not per-read sequences |
| UMI-tools | Only handles UMIs at fixed read positions (first/last N bases) |
| `samtools consensus` | Collapses reads into a consensus; no per-read output |

### Key features

- Control over insertions at the boundaries
- TSV and FASTA output, in wide or long layout
- Pipeable (both input and output) => easy integration with `samtools` and other tools

## Installation

### With conda

```bash
git clone https://github.com/fravadona/xumi.git
cd xumi
conda create -n xumi -c conda-forge -c bioconda python=3.10 pysam
conda activate xumi
pip install --no-deps .
```

### Requirements

- Python ≥ 3.10 (tested on 3.10, compatible with 3.7+)
- [pysam](https://github.com/pysam-developers/pysam) ≥ 0.20

## Quick start

```bash
# Single region
xumi.py -r chr1:100-200 -O fasta mapped.bam > out.fasta

# Multiple regions
xumi.py -r chr1:100-200,chr1:1000-1100 mapped.sam > out.tsv

# Regions from a BED file
xumi.py -R regions.bed -O tsv-long mapped.cram

# Pipe from samtools
samtools view -u mapped.bam chr1 | xumi.py -r chr1:100-200
```

## Usage

```text
xumi [-h] [-V] [-H] [-r REGIONS] [-R in.bed]
     [-a] [-b {both,left,right,none}]
     [-O {tsv,tsv-long,fasta,fasta-long}]
     [-o FILE] [aln.sam|aln.bam|aln.cram]
```
> **Note:** At least one of `-r` or `-R` is required.

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

### Output formats

| Format | Layout | Description |
|---|---|---|
| `tsv` (default) | wide | One row per read, one column per region |
| `tsv-long` | long | One row per read × region pair |
| `fasta` | wide | One record per read, regions joined with `-` |
| `fasta-long` | long | One record per read × region pair |

> **Note:** Reads with no extracted sequence are omitted.

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
xumi -r chr1:10-19,chr1:20-29 --boundary-insertions right mapped.bam
```

| Value | Left boundary | Right boundary |
|---|---|---|
| `both` (default) | included | included |
| `left` | included | excluded |
| `right` | excluded | included |
| `none` | excluded | excluded |

## Examples

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

### Extract dual-UMIs at known reference positions with pre and post filtering

```bash
# UMIs are at chr1:101-120 and chr1:2001-2020 (expected sizes of around 20bp)

# pre-filter reads with samtools (adjust filters to your needs)
samtools view -u -h \
    -F 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY' \
    -e '1700 <= length(seq) && length(seq) <= 2300' \
    -q 30 \
    input.bam | # append 'chr1' to speed up if the BAM is sorted+indexed

# extract the UMI regions, including boundary insertions
xumi -H -r 'chr1:101-120,chr1:2001-2020' -b both -O tsv |

# keep only records where both UMIs are within expected size range
awk '
    {umi1_len = length($2); umi2_len = length($3)}
    17 <= umi1_len && umi1_len <= 23 && 17 <= umi2_len && umi2_len <= 23
' FS='\t' > umis.tsv
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

## License

[MIT](LICENSE)

## Author

Rafael NAVAZA  
Institut Pasteur, Plateforme de Cristallographie  
rnavaza@pasteur.fr

