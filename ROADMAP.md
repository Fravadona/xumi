# xumi Roadmap

## Planned Features

### FASTQ-long output format (`fastq-long`)

A new output format for per-read × region extraction with quality scores.

- One record per read × region pair, header format: `@qname::region`
- Quality string sliced in parallel with the sequence using the same `[q0, q1)`
  interval computed during extraction
- Reads with missing quality scores (SAM `*`) use Phred 1 (`"`) as placeholder,
  consistent with `samtools fastq -v 1` behavior
- Requires adding `query_qualities: str | None` to `ReadProjection`

### SAM flag emission (`--emit-flag`)

Opt-in option to append the SAM flag in hex to the read name in all output formats.

- Appends `::0xFFF` to `qname` in all output formats
- Users with unique qnames in input (e.g. after pre-filtering with
  `samtools view -F 0x900`) do not need this option

### Alignment annotation (`--mark-ins`, `--mark-del`)

Options to annotate insertions and deletions in the extracted sequence.
Both options are independent and can be used simultaneously.

**`--mark-ins [uc|lc]`**
- Bare `--mark-ins`: prefix each inserted base with `_`, preserving the
  inserted sequence — follows `samtools consensus --mark-ins` convention;
  when used with `fastq-long`, `_` is also inserted at the corresponding
  quality string position with quality `!` (Phred 0), keeping sequence and
  quality string lengths equal — non-standard FASTQ, documented as such
- `lc`: inserted bases in lowercase, aligned bases unchanged
- `uc`: inserted bases in uppercase, aligned bases unchanged
- Incompatible with `--aligned-only` (no insertions to mark)

> Note: `lc`/`uc` are inspired by `bcftools consensus` conventions.

**`--mark-del`**
- Represent each deleted base (`D`) and reference skip (`N`) within the region
  as `*`, instead of silently omitting them (e.g. a 10bp deletion emits
  `**********`)
- Enables distinguishing deletions and reference skips from absent coverage
  in the output; users working with spliced alignments should be aware that
  `D` and `N` are not distinguished
- When used with `fastq-long`, `*` markers are assigned quality `!` (Phred 0)
  in the quality string — non-standard FASTQ, documented as such

> Note: `*` follows SAM and `bcftools consensus --mark-del` convention.

**Implementation note:** both options reuse the existing block/CIGAR index
infrastructure (`block_slice_for_region`, `cigar_op_index` in `MatchBlock`).
Walking `blocks[i0:i1]` and building the annotated string segment by segment
is sufficient for both modes.

### Basic read filters

Commonly needed pre-filters to reduce the need for an upstream `samtools view`.

**`-F`, `--exclude-flags` `INT|HEX|OCT|NAME[,NAME...]`**
- Exclude reads with any of the specified FLAG bits set
- Accepts decimal, hex (`0x900`), octal (`0100`), or comma-separated flag
  names (`UNMAP,SECONDARY,SUPPLEMENTARY`) — case-insensitive, consistent
  with `samtools view -F`
- Supported flag names:
  `PAIRED`, `PROPER_PAIR`, `UNMAP`, `MUNMAP`, `REVERSE`, `MREVERSE`,
  `READ1`, `READ2`, `SECONDARY`, `QCFAIL`, `DUP`, `SUPPLEMENTARY`

**`-q`, `--min-MQ` `INT`**
- Exclude reads with mapping quality below `INT`
- Default: 0 (no filtering), consistent with `samtools view -q`

Both filters are applied before extraction, with no interaction with extraction
logic.

### Pure Python fallback (no pysam)

A fallback SAM parser for environments where pysam is unavailable (e.g. Windows,
restricted HPC environments).

- SAM-only (plain text); BAM/CRAM support requires pysam
- `ReadProjection.from_cigar()` is already pysam-independent by design; only
  the IO layer needs a fallback implementation
- `pyproject.toml`: pysam declared as optional on Windows via
  `pysam >= 0.20; sys_platform != 'win32'`
- No changes to extraction or output logic

## Internal Refactoring (non-breaking)

### `ReadProjection` factory methods

- `ReadProjection.from_cigar(query_sequence, cigar, reference_start, *, query_qualities, track_deletions)` —
  generic constructor accepting CIGAR string or tuples, SAM quality string or
  decoded list; lives in the generic section, no pysam dependency; planned,
  not yet implemented
- `ReadProjection.from_aln(aln)` — thin pysam adapter calling `from_cigar`;
  handles unmapped/missing SEQ/CIGAR guards; uses `aln.query_qualities_str`
  when available (pysam ≥ 0.22), falls back to manual encoding otherwise;
  planned, not yet implemented

### Code reorganization

Splitting the current "Application-Specific Functions" section into named
subsections for readability:

1. **Genomic Types** — `Region`, `parse_bed`, `parse_region`
2. **Sequence Utilities** — `FASTA_LINE_WIDTH`, `fold_sequence`
3. **CIGAR Projection** — `COPS`, `MatchBlock`, `ReadProjection`,
   `block_slice_for_region`
4. **Extraction** — `RegionExtractionMode`, insertion length helpers,
   `extract_aligned_bases`, `query_slice_for_region`, `extract_query_slice`,
   `extract_region_from_projection`
5. **Output** — `OutputLayout`, `OutputFormat`, `_OUTPUT_MODES`,
   `ExtractionResult`, `emit_wide`, `emit_long`
6. **CLI** — `Config`, `parse_args`
7. **Application** — `_emit_header`, `_extract_all`, `run`, `XUMI_DEBUG`, `main`
