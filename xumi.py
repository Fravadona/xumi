#!/usr/bin/env python3
from __future__ import annotations

# ############################################################################ #

"""
XUMI: Extract subsequences from read alignments at specified genomic regions,
      with full CIGAR awareness. Originally developed for inline UMI extraction.

Author:      Rafael NAVAZA <rnavaza@pasteur.fr>
Affiliation: Institut Pasteur, Plateforme de Cristallographie
License:     MIT (see LICENSE file)
URL:         https://github.com/fravadona/xumi
"""

# ############################################################################ #

__version__ = "1.0.2"

# ############################################################################ #

# --- Core libraries
import os
import re
import sys
import gzip
import signal
import argparse
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from collections import defaultdict
from collections.abc import Iterator
from bisect import bisect_left, bisect_right

# --- Third-party libraries
import pysam

# ############################################################################ #
#
#  Generic Library Functions
#
# ############################################################################ #

@dataclass(frozen=True)
class Region:
    """
    A genomic region.

    Attributes
    ----------
    chrom : str
        Chromosome name.
    start : int
        Start coordinate (0-based).
    stop  : int
        End coordinate (exclusive).

    Properties
    ----------
    __str__ : str
        A standard representation with 1-based coordinates incl.
        Example: 'chr1:101-200'
    length : int
        Region length in bases.
    """
    chrom: str
    start: int
    stop:  int

    @property
    def length(self):
        return self.stop - self.start

    def __str__(self):
        return f"{self.chrom}:{self.start+1}-{self.stop}"

# ---------------------------------------------------------------------------- #

def parse_bed(path: str | Path) -> Iterator[Region]:

    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open

    with opener(path, 'rt') as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\r\n")

            if not line or line.startswith(('#', 'track', 'browser')):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}: line {line_no}: "
                    f"invalid BED record format -- {repr(line)}"
                )

            # As BED is strictly tab-delimited, the first field could
            # theoretically contain significant whitespaces; in
            # practice it can't so using strip() is safe here.
            chrom = fields[0].strip()

            if not (chrom and fields[1].isascii() and fields[1].isdigit()
                          and fields[2].isascii() and fields[2].isdigit()):
                raise ValueError(
                    f"{path}: line {line_no}: "
                    f"illegal value(s) in BED record -- {repr(line)}"
                )

            start = int(fields[1])
            stop  = int(fields[2])

            if start < 0 or start >= stop:
                raise ValueError(
                    f"{path}: line {line_no}: "
                    f"invalid BED region boundaries -- {repr(line)}"
                )

            yield Region(chrom=chrom, start=start, stop=stop)

# ---------------------------------------------------------------------------- #

# SAMtools-style region string parsing.
#
# SAMtools documentation specifies the format as RNAME[:STARTPOS[-ENDPOS]]
# with 1-based coordinates. Here we restrict it to RNAME:START-END as it
# doesn't make sense to extract full contigs in a region extraction tool.
#
# Example:
#   chr1:100-200 -> chr1 from 100 to 200 (inclusive).
#

_REGION_REGEX = re.compile(r'(?P<rname>.+):(?P<start>\d+)-(?P<stop>\d+)$')

def parse_region(region_str: str) -> Region:
    m = _REGION_REGEX.match(region_str)
    if not m:
        raise ValueError(f"invalid region format: {repr(region_str)}")

    rname = m.group('rname')
    start = int(m.group('start')) - 1
    stop  = int(m.group('stop'))

    if start < 0 or start >= stop:
        raise ValueError(f"invalid region boundaries: {repr(region_str)}")

    return Region(chrom=rname, start=start, stop=stop)


# ############################################################################ #
#
#  Application-Specific Functions
#
# ############################################################################ #

# Alias for pysam.CIGAR_OPS
COPS = pysam.CIGAR_OPS

# ---------------------------------------------------------------------------- #

class RegionExtractionMode(Enum):
    ALIGNED_BASES_ONLY = auto()
    QUERY_SLICE        = auto()
    QUERY_SLICE_WITH_LEFT_BOUNDARY_INSERTIONS  = auto()
    QUERY_SLICE_WITH_RIGHT_BOUNDARY_INSERTIONS = auto()
    QUERY_SLICE_WITH_BOTH_BOUNDARY_INSERTIONS  = auto()

# ---------------------------------------------------------------------------- #

@dataclass(frozen=True) # we can add "slots=True" starting from Python >= 3.10
class MatchBlock:
    """
    One contiguous aligned block from a single CIGAR op of type M/=/X.

    Coordinates are 0-based, half-open.
    - Ref interval   : [ref_start, ref_stop)
    - Query interval : [query_start, query_stop)

    cigar_op_index is the index into alignment.cigartuples for this M/=/X op.
    """
    cigar_op_index: int
    ref_start: int
    ref_stop: int
    query_start: int
    query_stop: int

# ---------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ReadProjection:
    query_sequence: str
    cigar_tuples: list[tuple[int, int]]
    blocks: list[MatchBlock]
    block_ref_starts: list[int]
    block_ref_stops: list[int]

# ############################################################################ #

def make_read_projection(aln: pysam.AlignedSegment) -> ReadProjection | None:
    """
    Create a ReadProjection directly from a pysam.AlignedSegment.
    It computes the match blocks and bisect-indices in a single CIGAR pass.
    Returns None if the read can't be projected (unmapped or missing SEQ/CIGAR).
    """
    if aln.is_unmapped or aln.query_sequence is None or aln.cigartuples is None:
        return None

    qry_seq = aln.query_sequence
    cig_tup = aln.cigartuples

    qry_pos = 0
    ref_pos = aln.reference_start

    blocks: list[MatchBlock] = []
    block_ref_starts: list[int] = []
    block_ref_stops:  list[int] = []

    for op_idx, (op, op_len) in enumerate(cig_tup):
        if op in (COPS.CMATCH, COPS.CEQUAL, COPS.CDIFF):
            r0 = ref_pos
            r1 = ref_pos + op_len
            q0 = qry_pos
            q1 = qry_pos + op_len

            blocks.append(MatchBlock(
                cigar_op_index = op_idx,
                ref_start      = r0,
                ref_stop       = r1,
                query_start    = q0,
                query_stop     = q1,
            ))
            block_ref_starts.append(r0)
            block_ref_stops.append(r1)

            ref_pos = r1
            qry_pos = q1

        elif op in (COPS.CINS, COPS.CSOFT_CLIP):
            qry_pos += op_len

        elif op in (COPS.CDEL, COPS.CREF_SKIP):
            # create a zero-length query block to mark the deletion
            blocks.append(MatchBlock(
                cigar_op_index = op_idx,
                ref_start = ref_pos,
                ref_stop = ref_pos + op_len,
                query_start = qry_pos,
                query_stop = qry_pos
            ))
            block_ref_starts.append(ref_pos)
            block_ref_stops.append(ref_pos + op_len)

            ref_pos += op_len

        else:
            # H, P, etc.
            pass

    # invariant: CIGAR must consume entire query
    assert qry_pos == aln.query_length, (
        f"SEQ/CIGAR lengths mismatch: "
        f"{aln.query_name} -- {aln.query_length} ≠ {qry_pos}"
    )

    return ReadProjection(
        query_sequence   = qry_seq,
        cigar_tuples     = cig_tup,
        blocks           = blocks,
        block_ref_starts = block_ref_starts,
        block_ref_stops  = block_ref_stops,
    )

# ---------------------------------------------------------------------------- #

def block_slice_for_region(
    proj: ReadProjection,
    region: Region
) -> tuple[int, int] | None:
    """
    Returns the smallest [i0, i1) interval such as proj.blocks[i0:i1]
    overlaps the region completely; i1 is exclusive.
    Returns None if no block overlaps.
    """
    starts = proj.block_ref_starts
    stops  = proj.block_ref_stops

    i0 = bisect_right(stops, region.start)   # first with stop > region.start
    i1 = bisect_left(starts, region.stop)    # first with start >= region.stop

    return None if i0 >= i1 else (i0, i1)

# ############################################################################ #

def extract_aligned_bases(proj: ReadProjection, region: Region) -> str | None:
    """
    Return concatenation of aligned read bases (M/=/X) overlapping `region`.
    Excludes insertions. Returns None if there is no aligned overlap.
    """
    span = block_slice_for_region(proj, region)
    if span is None:
        return None

    i0, i1 = span
    blocks = proj.blocks
    qs = proj.query_sequence

    out: list[str] = []
    for i in range(i0, i1):
        b = blocks[i]
        # skip deletion/skip block
        if b.query_start == b.query_stop:
            continue
        ov_ref_start = max(b.ref_start, region.start)
        ov_ref_stop  = min(b.ref_stop,  region.stop)
        if ov_ref_start < ov_ref_stop:
            q0 = b.query_start + (ov_ref_start - b.ref_start)
            q1 = b.query_start + (ov_ref_stop  - b.ref_start)
            out.append(qs[q0:q1])

    # If blocks overlap the region, each must produce at least one base
    # (guaranteed by the bisect logic and the ov_ref_start < ov_ref_stop
    # guard), so `out` is non-empty whenever `span` is not None.
    return "".join(out) if out else None

# ############################################################################ #

def left_insertion_length(
    cig_tup: list[tuple[int, int]],
    cigar_op_index: int
) -> int:
    op_idx = cigar_op_index - 1
    n_ins = 0
    while op_idx >= 0 and cig_tup[op_idx][0] == COPS.CINS:
        n_ins += cig_tup[op_idx][1]
        op_idx -= 1
    return n_ins

# ---------------------------------------------------------------------------- #

def right_insertion_length(
    cig_tup: list[tuple[int, int]],
    cigar_op_index: int
) -> int:
    op_idx = cigar_op_index + 1
    n_ins = 0
    while op_idx < len(cig_tup) and cig_tup[op_idx][0] == COPS.CINS:
        n_ins += cig_tup[op_idx][1]
        op_idx += 1
    return n_ins

# ---------------------------------------------------------------------------- #

def query_slice_for_region(
    proj: ReadProjection,
    region: Region,
    *,
    include_left_boundary_insertions: bool,
    include_right_boundary_insertions: bool,
) -> tuple[int, int] | None:
    """
    Return [q0, q1) for slicing proj.query_sequence[q0:q1]
    in the QUERY_SLICE* extraction modes.

    Anchors are the first and last aligned bases (M/=/X)
    overlapping the region. Insertions between anchors are
    inherently included because the output is a contiguous
    read substring. Optionally, the slice is extended to
    include insertion ops (I) immediately adjacent to the
    boundary anchor(s).
    """
    span = block_slice_for_region(proj, region)
    if span is None:
        return None

    i0, i1 = span
    blocks = proj.blocks

    first_block = blocks[i0]
    last_block  = blocks[i1 - 1]

    # ----- Inner anchors -----
    left_overlap_ref_start = max(first_block.ref_start, region.start)
    right_overlap_ref_stop = min(last_block.ref_stop,   region.stop)

    #qry_start_inner = first_block.query_start + (left_overlap_ref_start - first_block.ref_start)
    qry_start_inner = qry_start_inner = max(first_block.query_start, first_block.query_start + (region.start - first_block.ref_start))
    qry_stop_inner  = min(last_block.query_stop, last_block.query_start + (region.stop - last_block.ref_start)) # exclusive

    qry_start = qry_start_inner
    qry_stop  = qry_stop_inner

    if qry_start == qry_stop:
        # region covered entirely by deletions or skips
        return None

    # ----- Boundary insertions -----
    if include_left_boundary_insertions:
        # Only include insertions to the left of a match block if the region's
        # anchor is at the very beginning of that block. Otherwise, the
        # insertion is not truly 'adjacent' to the region boundary.
        if left_overlap_ref_start == first_block.ref_start:
            qry_start -= left_insertion_length(proj.cigar_tuples, first_block.cigar_op_index)

    if include_right_boundary_insertions:
        # include right-adjacent insertions only if right anchor is at end of
        # its match op
        if right_overlap_ref_stop == last_block.ref_stop:
            qry_stop += right_insertion_length(proj.cigar_tuples, last_block.cigar_op_index)

    # invariant: slice bounds must be well-ordered and in range
    #assert (0 <= qry_start <= qry_start_inner < qry_stop_inner <= qry_stop <= len(proj.query_sequence))

    return qry_start, qry_stop

# ---------------------------------------------------------------------------- #

def extract_query_slice(
    proj: ReadProjection,
    region: Region,
    *,
    include_left_boundary_insertions: bool,
    include_right_boundary_insertions: bool,
) -> str | None:

    interval = query_slice_for_region(
        proj,
        region,
        include_left_boundary_insertions = include_left_boundary_insertions,
        include_right_boundary_insertions = include_right_boundary_insertions,
    )
    if interval is None:
        return None
    q0, q1 = interval
    return proj.query_sequence[q0:q1]

# ---------------------------------------------------------------------------- #

def extract_region_from_projection(
    proj: ReadProjection,
    region: Region,
    mode: RegionExtractionMode,
) -> str | None:

    if mode == RegionExtractionMode.ALIGNED_BASES_ONLY:
        return extract_aligned_bases(proj, region)

    if mode == RegionExtractionMode.QUERY_SLICE:
        return extract_query_slice(
            proj, region,
            include_left_boundary_insertions=False,
            include_right_boundary_insertions=False,
        )

    if mode == RegionExtractionMode.QUERY_SLICE_WITH_LEFT_BOUNDARY_INSERTIONS:
        return extract_query_slice(
            proj, region,
            include_left_boundary_insertions=True,
            include_right_boundary_insertions=False,
        )

    if mode == RegionExtractionMode.QUERY_SLICE_WITH_RIGHT_BOUNDARY_INSERTIONS:
        return extract_query_slice(
            proj, region,
            include_left_boundary_insertions=False,
            include_right_boundary_insertions=True,
        )

    if mode == RegionExtractionMode.QUERY_SLICE_WITH_BOTH_BOUNDARY_INSERTIONS:
        return extract_query_slice(
            proj, region,
            include_left_boundary_insertions=True,
            include_right_boundary_insertions=True,
        )

    raise RuntimeError(f"unhandled extraction mode: {mode!r}")

# ############################################################################ #

class OutputLayout(Enum):
    WIDE = auto()
    LONG = auto()

class OutputFormat(Enum):
    TSV   = auto()
    FASTA = auto()

_OUTPUT_MODES = {
    'tsv':        (OutputFormat.TSV,   OutputLayout.WIDE),
    'tsv-long':   (OutputFormat.TSV,   OutputLayout.LONG),
    'fasta':      (OutputFormat.FASTA, OutputLayout.WIDE),
    'fasta-long': (OutputFormat.FASTA, OutputLayout.LONG),
}

# ---------------------------------------------------------------------------- #

FASTA_LINE_WIDTH = 80

def fold_sequence(seq: str, width: int = FASTA_LINE_WIDTH) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

# ---------------------------------------------------------------------------- #

def emit_wide(out, fmt, qname, extractions, tgt_regions):
    if fmt == OutputFormat.TSV:
        cols = [''] * len(tgt_regions)
        for rgn_idx, rgn, seq in extractions:
            cols[rgn_idx] = seq
        print(f"{qname}\t" + "\t".join(cols), file=out)

    elif fmt == OutputFormat.FASTA:
        seqs = [''] * len(tgt_regions)
        for rgn_idx, rgn, seq in extractions:
            seqs[rgn_idx] = seq
        print(f">{qname}", file=out)
        # The best compromise between "standard" FASTA and a meaningful output
        # is to join the sequences of the regions with a hyphen because they
        # they can't appear in SAM sequences and the user will be able
        # to retrieve them even when using a standard FASTA parser.
        print(fold_sequence('-'.join(seqs)), file=out)

# ---------------------------------------------------------------------------- #

def emit_long(out, fmt, qname, extractions):
    for rgn_idx, rgn, seq in extractions:
        if fmt == OutputFormat.TSV:
            print(f"{qname}\t{rgn}\t{seq}", file=out)

        elif fmt == OutputFormat.FASTA:
            print(f">{qname}::{rgn}", file=out)
            print(fold_sequence(seq), file=out)


# ############################################################################ #

# choices for CLI options; first-in-tuple is the default (when applicable)
_BOUNDARY_INSERTION_CHOICES = ("both", "left", "right", "none")
_OUTPUT_FORMAT_CHOICES = tuple(_OUTPUT_MODES)

# ---------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Config:
    tgt_regions: list[Region]
    extract_mode: RegionExtractionMode
    out_fmt: OutputFormat
    out_layout: OutputLayout
    no_header: bool
    output_path: str
    alignments_file: str

# ---------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> Config:
    cli = argparse.ArgumentParser(
        description = (
            "extract per‑read subsequences corresponding to "
            "specified reference regions."
        ),
        epilog = (
            "examples:\n"
            "  xumi.py -r chr1:100-200 mapped.bam\n"
            "  xumi.py -R regions.bed --aligned-only mapped.bam\n"
            "  xumi.py -R regions.bed --boundary-insertions left mapped.bam\n"
            "  samtools view -u mapped.bam chr1 | xumi.py -r chr1:100-200,chr1:1000-1100"
        ),
        formatter_class = argparse.RawTextHelpFormatter,
    )
    cli.add_argument(
        '-V', '--version',
        action = 'version',
        version = f'%(prog)s {__version__}',
    )
    cli.add_argument(
        '-H', '--no-header',
        action = 'store_true',
        default = False,
        help = (
            "suppress header line(s) in output."
        ),
    )
    cli.add_argument(
        '-r', '--regions',
        metavar = 'REGIONS',
        default = None,
        help = (
            "comma-separated list of regions, each specified as\n"
            "RNAME:START-END with 1-based coordinates."
        ),
    )
    cli.add_argument(
        '-R', '--regions-file',
        metavar = 'in.bed',
        default = None,
        help = (
            "BED file containing the regions to extract."
        ),
    )
    cli.add_argument(
        '-a', '--aligned-only',
        action = 'store_true',
        help=(
            "restrict extraction to reference-aligned bases (M/=/X)\n"
            "within the target region(s). by default, the entire\n"
            "portions of reads spanning those regions are returned.\n"
            "incompatible with --boundary-insertions."
        ),
    )
    cli.add_argument(
        '-b', '--boundary-insertions',
        choices = _BOUNDARY_INSERTION_CHOICES,
        default = None, # the default is _BOUNDARY_INSERTION_CHOICES[0] but
                        # we need to know when the user provides this option
        help = (
            "include insertions flanking the region boundaries on\n"
            "the specified side. default: both.\n"
            "only contiguous (I) operations directly adjacent to the\n"
            "first-last aligned bases in the region(s) are included.\n"
            "for contiguous regions extraction, use 'left' or 'right'\n"
            "to avoid including boundary insertions in both adjacent\n"
            "regions."
        ),
    )
    cli.add_argument(
        '-O', '--output-format',
         choices = _OUTPUT_FORMAT_CHOICES,
         default = _OUTPUT_FORMAT_CHOICES[0],
         help=(
            "output format. default: %(default)s.\n"
            "\n"
            "  tsv         one line per read, one column per region.\n"
            "  tsv-long    one line per read × region pair.\n"
            "  fasta       one record per read, regions delimited\n"
            "              with '-' (gap character) on the sequence\n"
            "              line.\n"
            "  fasta-long  one record per read × region pair.\n"
            "              header format: >QNAME::REGION\n"
            "\n"
            "reads with no extracted sequence are omitted.\n"
            "in wide formats (tsv, fasta), missing regions are\n"
            "represented as empty fields (tsv) or empty strings\n"
            "between delimiters (fasta)."
        ),
    )
    cli.add_argument(
        '-o', '--output',
        metavar = 'FILE',
        default = '-',
        help = (
            "output file. default: stdout."
        ),
    )
    cli.add_argument(
        'alignments_file',
        metavar = "aln.sam|aln.bam|aln.cram",
        nargs = '?',
        default = "-",
        help = (
            "SAM/BAM/CRAM file containing mapped reads.\n"
            "default: stdin."
        ),
    )
    args = cli.parse_args()

    # ------------------------------------------------------------------------ #

    # Validate option combinations
    if args.boundary_insertions and args.aligned_only:
        cli.error("--boundary-insertions conflicts with --aligned-only")

    # Set boundary_insertions option to the "real" default if not set
    if args.boundary_insertions is None:
        args.boundary_insertions = _BOUNDARY_INSERTION_CHOICES[0]

    # Map CLI options to internal RegionExtractionMode
    if args.aligned_only:
        extract_mode = RegionExtractionMode.ALIGNED_BASES_ONLY
    else:
        _BI_MODE_MAP = {
            "none":  RegionExtractionMode.QUERY_SLICE,
            "left":  RegionExtractionMode.QUERY_SLICE_WITH_LEFT_BOUNDARY_INSERTIONS,
            "right": RegionExtractionMode.QUERY_SLICE_WITH_RIGHT_BOUNDARY_INSERTIONS,
            "both":  RegionExtractionMode.QUERY_SLICE_WITH_BOTH_BOUNDARY_INSERTIONS,
        }
        extract_mode = _BI_MODE_MAP[args.boundary_insertions]

    # Collect regions
    if args.regions is None and args.regions_file is None:
        cli.error("at least one of -r/--regions or -R/--regions-file is required.")

    tgt_regions: list[Region] = []
    if args.regions is not None:
        tgt_regions.extend(map(parse_region, args.regions.split(',')))
    if args.regions_file is not None:
        tgt_regions.extend(parse_bed(args.regions_file))

    # Collect output format
    out_fmt, out_layout = _OUTPUT_MODES[args.output_format]

    return Config(
        tgt_regions     = tgt_regions,
        extract_mode    = extract_mode,
        out_fmt         = out_fmt,
        out_layout      = out_layout,
        no_header       = args.no_header,
        output_path     = args.output,
        alignments_file = args.alignments_file,
    )

# ---------------------------------------------------------------------------- #

def run(cfg: Config) -> None:
    if not cfg.tgt_regions:
        raise ValueError("no region provided.")

    # NOTE: pysam handles SAM/BAM/CRAM inputs from stdin automatically
    with pysam.AlignmentFile(cfg.alignments_file, 'r') as sam_file:
        sam_header = sam_file.header.to_dict()
        ref_lengths = {sq['SN']: sq['LN'] for sq in sam_header.get('SQ', [])}

        # Validate all regions against the SAM header
        for rgn in cfg.tgt_regions:
            if rgn.chrom not in ref_lengths:
                raise ValueError(f"region '{rgn}' not found in SAM header.")
            if rgn.stop > ref_lengths[rgn.chrom]:
                raise ValueError(f"out-of-bounds region '{rgn}'.")

        # Index regions by chromosome for fast lookup during iteration
        regions_by_chrom: dict[str, list[tuple[int, Region]]] = defaultdict(list)
        for i, rgn in enumerate(cfg.tgt_regions):
            regions_by_chrom[rgn.chrom].append((i, rgn))

        # Define the output stream
        if cfg.output_path == '-':
            out = sys.stdout
        elif cfg.output_path.endswith('.gz'):
            out = gzip.open(cfg.output_path, 'wt', compresslevel=6)
        else:
            out = open(cfg.output_path, 'w')

        try:
            if not cfg.no_header:
                _emit_header(out, cfg)

            # NOTE: Supplementary and secondary alignments are processed
            # like any other mapped read. A qname may therefore appear multiple
            # times in the output.
            for aln in sam_file:
                extractions = _extract_all(aln, regions_by_chrom, cfg.extract_mode)

                # Skip empty output records
                if not extractions:
                    continue

                qname = aln.query_name or '*'
                if cfg.out_layout == OutputLayout.WIDE:
                    emit_wide(out, cfg.out_fmt, qname, extractions, cfg.tgt_regions)
                else:
                    emit_long(out, cfg.out_fmt, qname, extractions)
        finally:
            if out is not sys.stdout:
                out.close()

# ---------------------------------------------------------------------------- #

def _emit_header(out, cfg: Config):
    if cfg.out_fmt == OutputFormat.TSV:
        if cfg.out_layout == OutputLayout.WIDE:
            cols = ["qname"] + [str(r) for r in cfg.tgt_regions]
        else:
            cols = ["qname", "region", "sequence"]
        print("#" + "\t".join(cols), file=out)

    elif cfg.out_fmt == OutputFormat.FASTA:
        if cfg.out_layout == OutputLayout.WIDE:
            regions = [str(r) for r in cfg.tgt_regions]
            print(";regions: " + ",".join(regions), file=out)

# ---------------------------------------------------------------------------- #

def _extract_all(aln, regions_by_chrom, extract_mode):
    extractions = []
    if not aln.is_unmapped:
        ref_name = aln.reference_name
        if ref_name in regions_by_chrom:
            proj = make_read_projection(aln)
            if proj is not None:
                for rgn_idx, rgn in regions_by_chrom[ref_name]:
                    seq = extract_region_from_projection(proj, rgn, extract_mode)
                    if seq is not None:
                        extractions.append((rgn_idx, rgn, seq))
    return extractions


# ---------------------------------------------------------------------------- #

XUMI_DEBUG = os.environ.get("XUMI_DEBUG", "0") not in ("0", "")
XUMI_DEBUG = True

def main() -> int:
    # restore SIGPIPE behavior of UNIX-tools
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    try:
        cfg = parse_args()
        run(cfg)
        return 0

    except KeyboardInterrupt:
        return 130

    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        if XUMI_DEBUG:
            import traceback
            traceback.print_exc()
        return 1

# ############################################################################ #

if __name__ == '__main__':
    sys.exit(main())

