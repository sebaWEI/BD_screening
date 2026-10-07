from pathlib import Path

from Bio import SeqIO

from bssf.cli import _header_metadata
from bssf.models import Target
from bssf.pipeline import generate_windows
from bssf.resources import EXAMPLE_LETM1, EXAMPLE_NSD2, package_root


def test_letm1_minus_strand_windows_match_bed_header() -> None:
    fasta = package_root() / EXAMPLE_LETM1["path"]
    record = SeqIO.read(fasta, "fasta")
    header = _header_metadata(record.description)
    start, end = int(header["start"]), int(header["end"])
    assert header["strand"] == "-"
    assert header["chrom"] == EXAMPLE_LETM1["chrom"]
    assert start == EXAMPLE_LETM1["start"]
    assert end == EXAMPLE_LETM1["end"]
    assert len(record.seq) == end - start == 2945
    windows = generate_windows(
        Target(str(record.seq), record.id, header["chrom"], start, end, "-"),
        window_size=40,
    )
    first = windows.iloc[0]
    last = windows.iloc[-1]
    assert first["start"] == end - 40
    assert first["end"] == end
    assert first["target_sequence"] == str(record.seq[:40]).upper()
    assert last["start"] == start
    assert last["end"] == start + 40
    assert last["target_sequence"] == str(record.seq[-40:]).upper()


def test_nsd2_plus_strand_windows_match_bed_header() -> None:
    fasta = package_root() / EXAMPLE_NSD2["path"]
    record = SeqIO.read(fasta, "fasta")
    header = _header_metadata(record.description)
    start, end = int(header["start"]), int(header["end"])
    assert header["strand"] == "+"
    assert len(record.seq) == end - start
    windows = generate_windows(
        Target(str(record.seq), record.id, header["chrom"], start, end, "+"),
        window_size=40,
    )
    first = windows.iloc[0]
    last = windows.iloc[-1]
    assert first["start"] == start
    assert first["end"] == start + 40
    assert last["start"] == end - 40
    assert last["end"] == end
