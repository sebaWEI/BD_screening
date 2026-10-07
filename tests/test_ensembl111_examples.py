"""Compare bundled example FASTAs to Ensembl 111 genome sequence."""

from __future__ import annotations

from Bio import SeqIO
from Bio.Seq import Seq
import pytest

from bssf.cli import _header_metadata
from bssf.fetch import (
    EnsemblArchiveError,
    assert_ensembl_111,
    fetch_bed_sequence,
    fetch_gene_utr,
    fetch_region,
)
from bssf.models import Target
from bssf.pipeline import generate_windows
from bssf.resources import EXAMPLE_LETM1, EXAMPLE_NSD2, package_root


def _load_example(spec: dict) -> tuple[str, dict[str, str]]:
    record = SeqIO.read(package_root() / spec["path"], "fasta")
    header = _header_metadata(record.description)
    return str(record.seq).upper().replace("U", "T"), header


def _plus_genome(chrom: str, start: int, end: int) -> str:
    try:
        return fetch_region(chrom, start + 1, end, 1)
    except EnsemblArchiveError as exc:
        pytest.skip(str(exc))


def _minus_transcript(chrom: str, start: int, end: int) -> str:
    try:
        return fetch_bed_sequence(chrom, start, end, "-")
    except EnsemblArchiveError as exc:
        pytest.skip(str(exc))


@pytest.fixture(scope="module")
def ensembl_111() -> str:
    try:
        return assert_ensembl_111()
    except EnsemblArchiveError as exc:
        pytest.skip(str(exc))


@pytest.mark.integration
def test_letm1_matches_ensembl111_genome_reverse_complement(ensembl_111: str) -> None:
    fasta, header = _load_example(EXAMPLE_LETM1)
    start, end = int(header["start"]), int(header["end"])
    chrom = header["chrom"]
    assert header["strand"] == "-"

    plus = _plus_genome(chrom, start, end)
    minus = _minus_transcript(chrom, start, end)
    rc_plus = str(Seq(plus).reverse_complement())
    assert fasta == minus == rc_plus

    windows = generate_windows(
        Target(fasta, "LETM1", chrom, start, end, "-"), window_size=40
    )
    first, last = windows.iloc[0], windows.iloc[-1]
    first_plus = plus[int(first["start"]) - start : int(first["end"]) - start]
    last_plus = plus[int(last["start"]) - start : int(last["end"]) - start]
    assert first["target_sequence"] == str(Seq(first_plus).reverse_complement()) == fasta[:40]
    assert last["target_sequence"] == str(Seq(last_plus).reverse_complement()) == fasta[-40:]


@pytest.mark.integration
def test_nsd2_matches_ensembl111_plus_strand_genome(ensembl_111: str) -> None:
    fasta, header = _load_example(EXAMPLE_NSD2)
    start, end = int(header["start"]), int(header["end"])
    chrom = header["chrom"]
    assert header["strand"] == "+"

    plus = _plus_genome(chrom, start, end)
    assert fasta == plus

    windows = generate_windows(
        Target(fasta, "NSD2", chrom, start, end, "+"), window_size=40
    )
    first, last = windows.iloc[0], windows.iloc[-1]
    assert first["target_sequence"] == plus[:40]
    assert last["target_sequence"] == plus[-40:]


@pytest.mark.integration
def test_example_fastas_match_ensembl111_canonical_utr(ensembl_111: str) -> None:
    letm_fasta, _ = _load_example(EXAMPLE_LETM1)
    nsd_fasta, _ = _load_example(EXAMPLE_NSD2)
    try:
        letm = fetch_gene_utr("LETM1")
        nsd = fetch_gene_utr("NSD2")
    except EnsemblArchiveError as exc:
        pytest.skip(str(exc))
    assert letm.transcript_id == EXAMPLE_LETM1["transcript_id"]
    assert letm.chrom == EXAMPLE_LETM1["chrom"]
    assert letm.start == EXAMPLE_LETM1["start"]
    assert letm.end == EXAMPLE_LETM1["end"]
    assert letm.strand == "-"
    assert letm.sequence == letm_fasta
    assert nsd.transcript_id == EXAMPLE_NSD2["transcript_id"]
    assert nsd.chrom == EXAMPLE_NSD2["chrom"]
    assert nsd.start == EXAMPLE_NSD2["start"]
    assert nsd.end == EXAMPLE_NSD2["end"]
    assert nsd.strand == "+"
    assert nsd.sequence == nsd_fasta
    assert letm.annotation_release == ensembl_111
