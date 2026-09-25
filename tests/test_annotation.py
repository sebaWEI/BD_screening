from pathlib import Path

from bsst.annotation import build_utr_index, utr_regions, write_utr_index, load_utr_index


def _gff3(tmp_path: Path) -> Path:
    text = """##gff-version 3
chr1	HAVANA	exon	100	179	.	+	.	ID=e1;transcript_id=ENST00000000001.1;exon_number=1
chr1	HAVANA	exon	300	399	.	+	.	ID=e2;transcript_id=ENST00000000001.1;exon_number=2
chr1	HAVANA	five_prime_UTR	100	149	.	+	.	transcript_id=ENST00000000001.1
chr1	HAVANA	three_prime_UTR	350	399	.	+	.	transcript_id=ENST00000000001.1
chr1	HAVANA	exon	300	399	.	-	.	ID=m1;transcript_id=ENST00000000002.1;exon_number=1
chr1	HAVANA	exon	100	179	.	-	.	ID=m2;transcript_id=ENST00000000002.1;exon_number=2
chr1	HAVANA	five_prime_UTR	350	399	.	-	.	transcript_id=ENST00000000002.1
chr1	HAVANA	three_prime_UTR	100	149	.	-	.	transcript_id=ENST00000000002.1
"""
    path = tmp_path / "mini.gff3"
    path.write_text(text)
    return path


def test_plus_and_minus_utr_map_onto_spliced_mrna(tmp_path: Path) -> None:
    index = build_utr_index(_gff3(tmp_path))
    assert index["ENST00000000001.1"] == [(1, 50, "5'UTR"), (131, 180, "3'UTR")]
    assert index["ENST00000000002.1"] == [(1, 50, "5'UTR"), (131, 180, "3'UTR")]
    assert utr_regions(index, "ENST00000000001.1", 40, 60) == "5'UTR"
    assert utr_regions(index, "ENST00000000001.1", 51, 130) == ""
    assert utr_regions(index, "ENST00000000001.1", 120, 140) == "3'UTR"
    assert utr_regions(index, "ENST999.1", 1, 10) == ""


def test_utr_index_roundtrip(tmp_path: Path) -> None:
    destination = tmp_path / "utr.tsv"
    write_utr_index(_gff3(tmp_path), destination)
    loaded = load_utr_index(destination)
    assert utr_regions(loaded, "ENST00000000002.1", 1, 10) == "5'UTR"
