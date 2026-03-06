#!/usr/bin/env python3

import subprocess
from pathlib import Path

TEST_DATA = Path(__file__).parent / "test_data"
SAM_FILE  = str(TEST_DATA / "cigar_ops.sam")
XUMI      = str(Path(__file__).parent.parent / "xumi.py")
REGIONS   = "chr1:101-120,chr1:121-140"


def run_xumi(*extra_args):
    result = subprocess.run(
        ["python3", XUMI, "-r", REGIONS, *extra_args, SAM_FILE],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def expected(filename):
    return (TEST_DATA / filename).read_text()


class TestSnapshotTSV:

    def test_aligned_only(self):
        assert run_xumi("-O", "tsv", "-a") == expected("cigar_ops_aligned.tsv")

    def test_slice_both(self):
        assert run_xumi("-O", "tsv", "-b", "both") == expected("cigar_ops_slice_both.tsv")

    def test_slice_none(self):
        assert run_xumi("-O", "tsv", "-b", "none") == expected("cigar_ops_slice_none.tsv")

    def test_slice_left(self):
        assert run_xumi("-O", "tsv", "-b", "left") == expected("cigar_ops_slice_left.tsv")

    def test_slice_right(self):
        assert run_xumi("-O", "tsv", "-b", "right") == expected("cigar_ops_slice_right.tsv")
