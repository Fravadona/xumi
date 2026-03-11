#!/usr/bin/env python3
import sys
import pysam

CIGAR_MAP = {
    0: "M", 1: "I", 2: "D", 3: "N",
    4: "S", 5: "H", 6: "P", 7: "=",
    8: "X"
}

def parse_region(region):
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start - 1, end   # 0-based half-open


def visualize_read(read, regions):

    chrom = read.reference_name or "*"
    qname = read.query_name

    ref_row = []
    anchor_row = []
    read_row = []
    cigar_row = []
    region_rows = [[] for _ in range(len(regions))]

    pairs = read.get_aligned_pairs(matches_only=False, with_cigar=True)
    query_seq = read.query_sequence

    qpos_eff = 0 # trying to work-around pysam bug in qpos when there is P ops
    prev_cigar_char = None
    for qpos, rpos, cigop in pairs:

        # cigar row
        cigar_char = CIGAR_MAP.get(cigop, "?")

        # reference row
        if cigar_char in ("P", "H"):
            ref_char = " "
        elif cigar_char == "I":
            ref_char = "-"
        elif rpos is not None:
            ref_char = "A"
        else:
            ref_char = " "

        # anchor row
        if cigar_char in ("P", "H"):
            anchor_char = " "
        elif rpos is not None:
            anchor_char = "|"
        else:
            anchor_char = " "

        # read row
        if cigar_char in ("P", "H"):
            read_char = " "
        elif qpos is None:
            read_char = "-"
        elif qpos_eff < 0 or qpos_eff >= len(query_seq):
            read_char = "?"  # sanity check
        else:
            read_char = query_seq[qpos_eff]

        # region row
        region_chars = []
        for chrom, rstart, rstop in regions:
            if cigar_char in ("P", "H"):
                region_chars.append(" ")
            elif rpos is not None and rstart <= rpos < rstop:
                region_chars.append("^")
            else:
                region_chars.append(" ")

        # split row by OP type in output
        if prev_cigar_char != cigar_char and prev_cigar_char is not None:
            ref_row.append(" ")
            anchor_row.append(" ")
            read_row.append(" ")
            cigar_row.append(" ")
            for rgn_row in region_rows:
                rgn_row.append(" ")

        prev_cigar_char = cigar_char

        ref_row.append(ref_char)
        anchor_row.append(anchor_char)
        read_row.append(read_char)
        cigar_row.append(cigar_char)
        for i, rgn_row in enumerate(region_rows): 
            rgn_row.append(region_chars[i])

        # increment only if op consumes query
        if cigar_char in ("M", "=", "X", "I", "S"):
            qpos_eff += 1

    for i, rgn_row in enumerate(region_rows):
        if "^" in rgn_row:
            i = rgn_row.index("^")
            for j in range(len(rgn_row) - 1, -1, -1):
                if rgn_row[j] == "^":
                    break
            rgn_row[i:j+1] = ["^"] * (j - i + 1)

    print()
    print("".join(ref_row)   + f"    {chrom}")
    print("".join(anchor_row))
    print("".join(read_row)  + f"    {qname}")
    print("".join(cigar_row) +"    CIGAR")
    for i, rgn_row in enumerate(region_rows):
        chrom, rstart, rstop = regions[i]
        print("".join(rgn_row)+ f"    {chrom}:{rstart+1}-{rstop}")
    print()


def main():

    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <sam/bam/cram> chr:start-end[,chr:start-end...]")

    path = sys.argv[1]

    regions = [parse_region(rgn) for rgn in sys.argv[2].split(",")]

    with pysam.AlignmentFile(path, "r") as f:
        for read in f:
            visualize_read(read, regions)


if __name__ == "__main__":
    main()
