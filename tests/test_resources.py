from bssf.resources import (
    DBSNP_B157_GRCH38P14,
    ENSEMBL_REST,
    GENCODE_V45_TRANSCRIPTS,
    bundled_gencode_fasta,
    catalog,
    chrom_to_refseq,
    refseq_to_chrom,
)


def test_pinned_resource_identifiers() -> None:
    assert DBSNP_B157_GRCH38P14["assembly"] == "GRCh38.p14"
    assert DBSNP_B157_GRCH38P14["refseq_accession"] == "GCF_000001405.40"
    assert DBSNP_B157_GRCH38P14["url"].endswith("GCF_000001405.40.gz")
    assert "nlm.nih.gov" in DBSNP_B157_GRCH38P14["url"]
    assert DBSNP_B157_GRCH38P14["tbi_filename"] == "GCF_000001405.40.gz.tbi"
    assert DBSNP_B157_GRCH38P14["refseq_to_chrom"]["NC_000001.11"] == "1"
    assert DBSNP_B157_GRCH38P14["refseq_to_chrom"]["NC_012920.1"] == "MT"
    assert chrom_to_refseq("4") == "NC_000004.12"
    assert chrom_to_refseq("chr4") == "NC_000004.12"
    assert refseq_to_chrom("NC_000004.12") == "4"
    assert GENCODE_V45_TRANSCRIPTS["n_transcripts"] == 252930
    assert GENCODE_V45_TRANSCRIPTS["ensembl_version"] == "111"
    assert GENCODE_V45_TRANSCRIPTS["assembly"] == "GRCh38.p14"
    assert GENCODE_V45_TRANSCRIPTS["url"].endswith("release_45/gencode.v45.transcripts.fa.gz")
    assert GENCODE_V45_TRANSCRIPTS["archive_md5"] == "415f81cc2f111fd6dbff5723de0cfa0b"
    assert len(GENCODE_V45_TRANSCRIPTS["archive_sha256"]) == 64
    assert len(GENCODE_V45_TRANSCRIPTS["fasta_sha256"]) == 64
    assert ENSEMBL_REST["pinned"] is True
    assert ENSEMBL_REST["base_url"] == "https://e111.rest.ensembl.org"
    assert ENSEMBL_REST["ensembl_version"] == "111"
    ids = [item["id"] for item in catalog()]
    assert ids == [
        "dbsnp_b157_grch38p14_tabix",
        "gencode_v45_transcripts_chr",
        "gencode_v45_annotation_chr_gff3",
        "ensembl_rest_111",
        "blastn_short",
        "example_letm1",
        "example_nsd2",
    ]


def test_bundled_gencode_fasta_matches_v45_if_present() -> None:
    path = bundled_gencode_fasta()
    if not path.is_file():
        return
    first = ""
    n_seq = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                n_seq += 1
                if not first:
                    first = line[1:].strip()
    assert first == GENCODE_V45_TRANSCRIPTS["first_header"]
    assert n_seq == GENCODE_V45_TRANSCRIPTS["n_transcripts"]
