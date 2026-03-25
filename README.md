# xumi

![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![GitHub](https://img.shields.io/github/release/fravadona/xumi)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18905757.svg)](https://doi.org/10.5281/zenodo.18905757)

CIGAR-aware extraction of per-read subsequences at specified reference regions from SAM/BAM/CRAM files.

## Overview

Originally developed for inline dual-UMI extraction,
xumi generalizes to any region-based per-read subsequence extraction from aligned reads.

### Motivation

Existing tools address related problems but leave a gap for per-read, region-specific extraction:

| Tool | Limitation |
|---|---|
| `samtools view` | Extracts whole reads, not the subsequences within region(s) |
| `bedtools getfasta` | Extracts reference sequences, not per-read sequences |
| `UMI-tools` | Extracts UMIs from patterns before alignment |
| `samtools consensus` | Collapses reads into a consensus, loses per-read information |
| `cigarillo (R)` | R library only, not a command-line utility |

**xumi** fills this gap and offers the following capabilities:

- Extraction of multiple regions per read
- CIGAR-aware control over insertions at region boundaries
- Output in TSV or FASTA, in wide or long layout
- Pipeable via stdin/stdout for easy integration with `samtools` and other tools
- Efficient processing of millions of reads

## Installation

### Requirements

- Python ≥ 3.10
- [pysam](https://github.com/pysam-developers/pysam) ≥ 0.20

### With conda

```bash
conda create -n xumi -c conda-forge -c bioconda xumi
conda activate xumi
```

## Quick start

```bash
# Extract subsequence in FASTA-format (default is TSV)
xumi -r chr1:100-200 -O fasta mapped.bam
```

```text
;regions: chr1:100:200
>read1
ATGCCGTTAAGCGGCGTACGTGCCCGATAGAGAGCTTACGACAGTGTACACGATGCCCGATCTAGTCAACCGGACTACGA
ACGTACGTTGCATGCA
>read2
...
```


## Usage

```text
xumi [options] <aln.{sam|bam|cram}>
```

### Options

| Option | Description |
|---|---|
| `-r`, `--regions` | Comma-separated regions: `chr1:100-200` (100 inclusive, 200 inclusive) |
| `-R`, `--regions-file` | BED file containing regions |
| `-a`, `--aligned-only` | Only output reference-aligned bases, *ignoring insertions within the region(s)* |
| `-b`, `--boundary-insertions` | Include boundary insertions: `both` (default), `left`, `right`, `none` |
| `-O`, `--output-format` | Output format (see below) |
| `-o`, `--output` | Output file (default: stdout, `.gz` for gzipped) |
| `-H`, `--no-header` | Suppress header line |
| `-h`, `--help` | Show help message |
| `-V`, `--version` | Show version |
> **Note:** At least one of `-r` or `-R` must be specified.

### Output formats

| Format | Layout | Description |
|---|---|---|
| `tsv` (default) | wide | One row per read, one column per region |
| `tsv-long` | long | One row per read × region pair |
| `fasta` | wide | One record per read, regions joined with `-` |
| `fasta-long` | long | One record per read × region pair |

> **Note:** Reads with no extracted sequence are omitted in output.

#### Example outputs (with header):

**TSV (wide)**

```text
#qname   chr1:1-8   chr1:101-108
read1    ACGTACGT   TGCATGCA
read2    ACGTACG    TGCATGCAT
```

**TSV (long)**

```text
#qname   region        sequence
read1    chr1:1-8      ACGTACGT
read1    chr1:101-108  TGCATGCA
read2    chr1:1-8      ACGTACG
read2    chr1:101-108  TGCATGCAT
```

**FASTA (wide)**

```text
;regions: chr1:1-8,chr1:101-108
>read1
ACGTACGT-TGCATGCA
>read2
ACGTACG-TGCATGCAT
```

**FASTA (long)**

```text
>read1::chr1:1-8
ACGTACGT
>read1::chr1:101-108
TGCATGCA
>read2::chr1:1-8
ACGTACG
>read2::chr1:101-108
TGCATGCAT
```

### Extraction modes (worked example)

Consider the following alignment:

```text
                         1             7
Chr1:                 -- AAA -- AAA -- C --
                         |||    |||    |
Read:            GGGG TT CCC TT AAA TT - TT
CIGAR:           SSSS II MMM II MMM II D II
```

Suppose we extract the region `Chr1:1-7`, which corresponds to `AAAAAAC` in reference space.
Here's what **xumi** yields using different modes:

| Mode | Extracted sequence |
|---|---|
| Aligned-only: `--aligned-only`               | CCCAAA |
| None: `--boundary-insertions none`           |   CCCTTAAATT |
| Left: `--boundary-insertions left`           | TTCCCTTAAATT |
| Right: `--boundary-insertions right`         |   CCCTTAAATTTT |
| Both (default): `--boundary-insertions both` | TTCCCTTAAATTTT |

> **Note:** Soft-clipped bases (`S`) at region boundaries are ignored, and deletions (`D`) serve
> as anchors only when they fall inside the region and the region also contains at least one mapped base.

## Advanced example

### Extract dual-UMIs at known reference positions with pre and post filtering

```bash
# UMIs are at chr1:101-120 and chr1:2001-2020 (expected sizes of around 20bp)

# pre-filter reads with samtools (adjust filters to your needs)
samtools view -u -h \
    -F 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY' \
    -e '1700 <= length(seq) && length(seq) <= 2300' \
    -q 30 \
    input.bam | # restrict to 'chr1' to speed up processing if BAM is sorted and indexed

# extract the UMI regions, including boundary insertions
xumi -H -r 'chr1:101-120,chr1:2001-2020' -b both -O tsv |

# post-filter sequences: only keep record if both UMIs are within size range (17-23 bp)
awk '
    {umi1_len = length($2); umi2_len = length($3)}
    17 <= umi1_len && umi1_len <= 23 && 17 <= umi2_len && umi2_len <= 23
' FS='\t' > extracted_umis.tsv
```

## Feedback and Issues

Your feedback, bug reports, and feature suggestions are highly welcome and help improve `xumi`.
Please use the [GitHub Issues page](https://github.com/fravadona/xumi/issues) for all contributions.

When reporting an issue, please include:
*   A clear description of the problem or proposed feature.
*   The `xumi` command(s) and input data (or minimal example) that caused the issue.
*   The output of `xumi -V` and your `pysam` version.
*   **For bugs:** Enable detailed logging with `XUMI_DEBUG=1 xumi ...` to provide a stack trace if an error occurs.

## License

[MIT](LICENSE)

## Author

Rafael NAVAZA (Institut Pasteur)

## Citation

If you use **xumi** in your work, please cite:

Navaza R. (2026). xumi: CIGAR-aware subsequence extraction from aligned reads. Zenodo. [DOI: 10.5281/zenodo.18905757](https://doi.org/10.5281/zenodo.18905757)

