import hashlib
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from bssf.models import Target
from bssf.db import download_file
from bssf.pipeline import (
    blast_hit_is_risk,
    blast_hit_is_self,
    complexity_failure,
    generate_windows,
    read_variant_vcf,
    variant_overlaps,
)


def test_windows_positive_and_negative() -> None:
    positive = generate_windows(Target("AACCGG", "x", "chr1", 10, 16, "+"), 3)
    negative = generate_windows(Target("AACCGG", "x", "1", 10, 16, "-"), 3)
    assert positive.loc[0, ["start", "end"]].tolist() == [10, 13]
    assert positive.loc[0, "bd_sequence"] == "GTT"
    assert negative.loc[0, ["start", "end"]].tolist() == [13, 16]
    assert negative.loc[1, ["start", "end"]].tolist() == [12, 15]
    assert generate_windows(Target("ACG", "short"), 40).empty


def test_complexity() -> None:
    assert complexity_failure("AAAAACGT") == "homopolymer"
    assert complexity_failure("ATATATAT") == "gc_out_of_range"
    assert complexity_failure("ACGTACGT") is None


def test_vcf_coordinates_chrom_and_af(tmp_path: Path) -> None:
    vcf = tmp_path / "small.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t11\t.\tA\tG\t.\t.\tAF=0.05\n"
        "1\t20\t.\tAT\tA\t.\t.\tAF=0.001\n"
        "4\t30\t.\tC\tT\t.\t.\tCAF=0.90,0.10;COMMON=1\n"
        "5\t40\t.\tG\tA\t.\t.\tFREQ=KOREAN:0.99,0.01|dbGaP_PopFreq:1,0;COMMON\n"
    )
    variants = read_variant_vcf(vcf, min_af=0.01)
    assert variants == [
        ("1", 10, 11, 0.05),
        ("4", 29, 30, 0.10),
        ("5", 39, 40, 0.01),
    ]
    assert variant_overlaps("1", 10, 11, variants)
    assert not variant_overlaps("chr1", 11, 12, variants)


def test_parse_refseq_common_via_local_scan(tmp_path: Path) -> None:
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "NC_000001.11\t100\trs1\tA\tG\t.\t.\tRS=1;COMMON\n"
        "NC_000001.11\t200\trs2\tC\tT\t.\t.\tRS=2\n"
        "NC_012920.1\t400\trs4\tA\tC\t.\t.\tRS=4;FREQ=x:0.9,0.1;COMMON\n"
    )
    all_sites = read_variant_vcf(vcf, common_only=False)
    assert ("1", 99, 100, None) in all_sites
    common = read_variant_vcf(vcf, common_only=True)
    assert common == [("1", 99, 100, None), ("mt", 399, 400, 0.1)]


def test_vcf_region_keeps_only_overlapping_chrom(tmp_path: Path) -> None:
    vcf = tmp_path / "small.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t11\t.\tA\tG\t.\t.\tAF=0.05\n"
        "4\t30\t.\tC\tT\t.\t.\tCAF=0.90,0.10;COMMON=1\n"
        "4\t1000\t.\tA\tG\t.\t.\tCAF=0.50,0.50;COMMON=1\n"
    )
    variants = read_variant_vcf(vcf, chrom="4", start=20, end=40)
    assert variants == [("4", 29, 30, 0.10)]


def test_tabix_region_query(tmp_path: Path) -> None:
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        return
    vcf = tmp_path / "region.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "4\t30\t.\tC\tT\t.\t.\tCAF=0.90,0.10;COMMON=1\n"
        "4\t500\t.\tA\tG\t.\t.\tCAF=0.50,0.50;COMMON=1\n"
    )
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    gz = tmp_path / "region.vcf.gz"
    subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
    variants = read_variant_vcf(gz, chrom="4", start=20, end=40)
    assert variants == [("4", 29, 30, 0.10)]


def test_blast_boundary_identity_coverage_and_exact_self() -> None:
    hit = pd.Series(
        {"sseqid": "ref|NM_001.2|", "sallacc": "NM_001.2", "stitle": "GENE1 transcript",
         "length": 20, "pident": 90}
    )
    short = pd.Series(
        {"sseqid": "ref|NM_002.1|", "sallacc": "NM_002.1", "stitle": "GENE2 transcript",
         "length": 19, "pident": 100}
    )
    assert blast_hit_is_risk(hit, query_length=20)
    assert blast_hit_is_risk(hit, query_length=40)
    assert not blast_hit_is_risk(short, query_length=40)
    assert not blast_hit_is_risk(
        pd.Series({"sseqid": "x", "sallacc": "", "stitle": "GENE2", "length": 15, "pident": 100}),
        query_length=40,
    )
    assert not blast_hit_is_risk(
        hit, query_length=40, min_identity=80.0, min_coverage=0.8
    )
    assert blast_hit_is_self(hit, ["NM_001.2"])
    assert blast_hit_is_self(hit, ["NM_001"])
    assert blast_hit_is_self(hit, ["GENE1"])
    assert not blast_hit_is_self(hit, ["GENE"])
    gencode = pd.Series(
        {
            "sseqid": (
                "ENST00000302787.3|ENSG00000168924.15|-|"
                "LETM1-201|LETM1|5371|protein_coding|"
            ),
            "sallacc": "",
            "stitle": "",
            "length": 40,
            "pident": 100,
        }
    )
    assert blast_hit_is_self(gencode, ["ENST00000302787"])
    assert blast_hit_is_self(gencode, ["LETM1"])
    assert not blast_hit_is_self(gencode, ["NSD2"])


def test_stream_download_and_checksum(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_bytes(b"##fileformat=VCFv4.2\n")
    checksum = hashlib.sha256(source.read_bytes()).hexdigest()
    destination = tmp_path / "db" / "variants.vcf"
    assert download_file(
        source.as_uri(), destination, expected_sha256=checksum
    ).read_bytes() == source.read_bytes()


def test_download_rejects_checksum_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_bytes(b"##fileformat=VCFv4.2\n")
    destination = tmp_path / "out.vcf"
    try:
        download_file(source.as_uri(), destination, expected_sha256="0" * 64)
    except ValueError as exc:
        assert "checksum mismatch" in str(exc)
    else:
        raise AssertionError("expected checksum mismatch")
