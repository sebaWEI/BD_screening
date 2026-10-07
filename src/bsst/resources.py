"""Pinned scientific resources used by bsst.

Names, URLs, and version strings are the identifiers assigned by the data
producer. Do not rename a file's contents after another database.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# NCBI dbSNP
# Header inside the file: ##source=dbSNP ##dbSNP_BUILD_ID=151 ##reference=GRCh38.p7
# ---------------------------------------------------------------------------
DBSNP_B151_GRCH38P7_COMMON_ALL: dict[str, Any] = {
    "id": "dbsnp_b151_grch38p7_common_all",
    "role": "variant_vcf",
    "producer": "NCBI dbSNP",
    "filename": "dbSNP_b151_GRCh38p7_common_all_20180418.vcf.gz",
    "url": (
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/"
        "VCF/common_all_20180418.vcf.gz"
    ),
    "directory_url": (
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/"
    ),
    "assembly": "GRCh38.p7",
    "variant_source": "NCBI dbSNP",
    "variant_release": "b151 fileDate=20180418 common_all",
    "chrom_style": "unprefixed",  # 1, 2, ..., X, Y, MT — not chr1
    "frequency_fields": ("CAF", "TOPMED"),
    "header_checks": (
        "##source=dbSNP",
        "##dbSNP_BUILD_ID=151",
        "##reference=GRCh38.p7",
        "##fileDate=20180418",
    ),
    # SHA-256 of common_all_20180418.vcf.gz. NCBI publishes MD5 only;
    # that digest matches this file: a274dcecff9cfe6084eaef848080ad8d
    "sha256": "a4b5ca65be831ef071891c8f6eb4925a83a1a9ecaa84a2585b0b71c7b217d564",
    "md5": "a274dcecff9cfe6084eaef848080ad8d",
    "md5_url": (
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/"
        "VCF/common_all_20180418.vcf.gz.md5"
    ),
}

# ---------------------------------------------------------------------------
# GENCODE transcript sequences (BLAST subject)
# Identified from the local FASTA: 252930 records; protein_coding=89110;
# nonsense_mediated_decay=21427; lncRNA=57722; first header DDX11L2-202.
# Those counts match GENCODE Release 45 CHR statistics exactly.
# ---------------------------------------------------------------------------
GENCODE_V45_TRANSCRIPTS: dict[str, Any] = {
    "id": "gencode_v45_transcripts_chr",
    "role": "blast_subject",
    "producer": "GENCODE",
    "filename": "gencode.v45.transcripts.fa",
    "archive_filename": "gencode.v45.transcripts.fa.gz",
    "url": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_45/gencode.v45.transcripts.fa.gz"
    ),
    "directory_url": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/"
    ),
    "release_page": "https://www.gencodegenes.org/human/release_45.html",
    "stats_page": "https://www.gencodegenes.org/human/stats_45.html",
    "transcriptome_source": "GENCODE",
    "transcriptome_release": "GENCODE 45 / Ensembl 111 / 2024-01 / CHR transcripts",
    "release": "45",
    "ensembl_version": "111",
    "assembly": "GRCh38.p14",
    "regions": "CHR",  # reference chromosomes including MT; not patches/haplotypes
    "n_transcripts": 252930,
    "n_protein_coding_transcripts": 89110,
    "n_nmd_transcripts": 21427,
    "n_lncrna_transcripts": 57722,
    "blast_db_name": "gencode_v45_transcripts_db",
    "blast_title": "GENCODE v45 transcripts CHR GRCh38.p14",
    "header_format": (
        "transcript_id|gene_id|havana_gene_id|havana_transcript_id|"
        "transcript_name|gene_name|length|transcript_biotype|"
    ),
    "first_header": (
        "ENST00000456328.2|ENSG00000290825.1|-|OTTHUMT00000362751.1|"
        "DDX11L2-202|DDX11L2|1657|lncRNA|"
    ),
    "last_header": (
        "ENST00000387461.2|ENSG00000210196.2|-|-|MT-TP-201|MT-TP|68|Mt_tRNA|"
    ),
    # SHA-256 of gencode.v45.transcripts.fa.gz. GENCODE MD5SUMS:
    # 415f81cc2f111fd6dbff5723de0cfa0b  gencode.v45.transcripts.fa.gz
    "archive_sha256": (
        "43bc5f0eab276b3a765a42186d92a4318138a68fd2f60ee9fd2dd1104d02bde2"
    ),
    "archive_md5": "415f81cc2f111fd6dbff5723de0cfa0b",
    "md5sums_url": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_45/MD5SUMS"
    ),
    # SHA-256 of the uncompressed FASTA (gzip -dc of the archive).
    "fasta_sha256": (
        "4263185345dbc8da8f3092b314671b6d0b0bea79ba7fe4b58794f609dd78c998"
    ),
}

# CHR GFF3 from the same release. five_prime_UTR / three_prime_UTR are genomic;
# db init projects them onto spliced transcript coordinates.
GENCODE_V45_ANNOTATION: dict[str, Any] = {
    "id": "gencode_v45_annotation_chr_gff3",
    "role": "utr_annotation",
    "producer": "GENCODE",
    "filename": "gencode.v45.annotation.gff3.gz",
    "utr_index_filename": "gencode.v45.utr_on_transcript.tsv",
    "url": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_45/gencode.v45.annotation.gff3.gz"
    ),
    "release": "45",
    "ensembl_version": "111",
    "assembly": "GRCh38.p14",
    "transcriptome_release": "GENCODE 45 CHR GFF3 → spliced 5′UTR / 3′UTR",
    "archive_md5": "e17bf2c2d47a0cdf28f62591fb4600ed",
    "archive_sha256": "",
}

# ---------------------------------------------------------------------------
# Ensembl REST archive matching GENCODE 45 (Ensembl 111, January 2024)
# https://e111.rest.ensembl.org is the versioned archive, not rest.ensembl.org.
# Archives are typically kept ~5 years; FASTA remains the durable interface.
# ---------------------------------------------------------------------------
ENSEMBL_REST: dict[str, Any] = {
    "id": "ensembl_rest_111",
    "role": "canonical_3utr_fetch",
    "producer": "EMBL-EBI Ensembl",
    "base_url": "https://e111.rest.ensembl.org",
    "docs_url": "https://e111.rest.ensembl.org",
    "live_url": "https://rest.ensembl.org",
    "default_species": "homo_sapiens",
    "ensembl_version": "111",
    "release": "111",
    "gencode_release": "45",
    "assembly": "GRCh38.p14",
    "pinned": True,
    "endpoints": {
        "lookup_symbol": "/lookup/symbol/{species}/{gene}?expand=1",
        "sequence_region": "/sequence/region/{species}/{chrom}:{start}..{end}:{strand}",
        "info_software": "/info/software",
        "info_data": "/info/data",
    },
    "note": (
        "Convenience fetcher pinned to the Ensembl 111 REST archive so it "
        "matches GENCODE 45. Pass --utr with a stored 3'UTR for a durable run. "
        "If this archive is retired, use a stored FASTA; do not silently "
        "fall back to rest.ensembl.org."
    ),
}

# ---------------------------------------------------------------------------
# Native tools (not Python packages; not installed by uv)
# ---------------------------------------------------------------------------
BLASTN_SHORT: dict[str, Any] = {
    "id": "blastn_short",
    "role": "offtarget_search",
    "producer": "NCBI BLAST+",
    "program": "blastn",
    "docs_url": "https://www.ncbi.nlm.nih.gov/books/NBK279690/",
    "task": "blastn-short",
    "dbtype": "nucl",
    "max_target_seqs": 100,
    "query_alphabet": "DNA",  # T, never U
    "default_min_length": 20,
}

# ---------------------------------------------------------------------------
# Bundled example 3′UTRs (coordinates are 0-based BED)
# ---------------------------------------------------------------------------
EXAMPLE_LETM1: dict[str, Any] = {
    "id": "example_letm1",
    "role": "example_target",
    "producer": "bsst examples (Ensembl-style 3′UTR FASTA)",
    "path": "examples/LETM1.fasta",
    "gene": "LETM1",
    "transcript_id": "ENST00000302787",
    "gencode_v45_record": (
        "ENST00000302787.3|ENSG00000168924.15|OTTHUMG00000121149.5|"
        "OTTHUMT00000241634.2|LETM1-201|LETM1|5371|protein_coding|"
    ),
    "chrom": "4",
    "start": 1811478,
    "end": 1814423,
    "strand": "-",
    "verified_against": "Ensembl REST archive 111 sequence/region and lookup/symbol",
}

EXAMPLE_NSD2: dict[str, Any] = {
    "id": "example_nsd2",
    "role": "example_target",
    "producer": "bsst examples (Ensembl-style 3′UTR FASTA)",
    "path": "examples/NSD2.fasta",
    "gene": "NSD2",
    "transcript_id": "ENST00000508803",
    "gencode_v45_record": (
        "ENST00000508803.6|ENSG00000109685.19|OTTHUMG00000121147.15|"
        "OTTHUMT00000366357.3|NSD2-218|NSD2|7560|protein_coding|"
    ),
    "chrom": "4",
    "start": 1978909,
    "end": 1982192,
    "strand": "+",
    "verified_against": "Ensembl REST archive 111 sequence/region and lookup/symbol",
}


def catalog() -> tuple[dict[str, Any], ...]:
    return (
        DBSNP_B151_GRCH38P7_COMMON_ALL,
        GENCODE_V45_TRANSCRIPTS,
        GENCODE_V45_ANNOTATION,
        ENSEMBL_REST,
        BLASTN_SHORT,
        EXAMPLE_LETM1,
        EXAMPLE_NSD2,
    )


def package_root() -> Path:
    """Repository root in an editable src/ layout; otherwise CWD."""
    here = Path(__file__).resolve()
    repo = here.parents[2]
    if (repo / "pyproject.toml").is_file() and (repo / "src" / "bsst").is_dir():
        return repo
    return Path.cwd()


def data_dir() -> Path:
    return package_root() / "data"


def bundled_dbsnp_vcf() -> Path:
    return data_dir() / DBSNP_B151_GRCH38P7_COMMON_ALL["filename"]


def bundled_gencode_fasta() -> Path:
    return data_dir() / GENCODE_V45_TRANSCRIPTS["filename"]


def bundled_gencode_gff3() -> Path:
    return data_dir() / GENCODE_V45_ANNOTATION["filename"]


def bundled_utr_index() -> Path:
    return data_dir() / GENCODE_V45_ANNOTATION["utr_index_filename"]


def discover_utr_index() -> Path | None:
    path = bundled_utr_index()
    return path if path.is_file() else None


def bundled_blast_db_prefix() -> Path:
    return data_dir() / GENCODE_V45_TRANSCRIPTS["blast_db_name"]


def blast_db_is_present(prefix: str | Path) -> bool:
    path = Path(prefix)
    return path.with_suffix(".nhr").is_file() or path.with_suffix(".nin").is_file()


def discover_variant_vcf() -> Path | None:
    path = bundled_dbsnp_vcf()
    return path if path.is_file() else None


def discover_blast_db() -> str | None:
    prefix = bundled_blast_db_prefix()
    return str(prefix) if blast_db_is_present(prefix) else None


def tabix_index_path(vcf: Path) -> Path:
    return Path(str(vcf) + ".tbi")


def verify_dbsnp_vcf(path: Path) -> tuple[bool, str]:
    """Check producer header fields inside the VCF, not the file name."""
    opener = gzip.open if path.name.endswith(".gz") else open
    headers: list[str] = []
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                headers.append(line.rstrip("\n"))
                if len(headers) >= 40:
                    break
    except OSError as exc:
        return False, str(exc)
    text = "\n".join(headers)
    missing = [
        field for field in DBSNP_B151_GRCH38P7_COMMON_ALL["header_checks"]
        if field not in text
    ]
    if missing:
        return False, "missing " + ", ".join(missing)
    return True, "dbSNP b151 GRCh38.p7 common_all"


def verify_gencode_fasta(path: Path) -> tuple[bool, str]:
    try:
        with path.open(encoding="utf-8") as handle:
            first = handle.readline().rstrip("\n")
    except OSError as exc:
        return False, str(exc)
    expected = ">" + GENCODE_V45_TRANSCRIPTS["first_header"]
    if first != expected:
        return False, first[:120] or "empty FASTA"
    return True, "GENCODE 45 CHR first header"


def verify_blast_db(prefix: str | Path) -> tuple[bool, str]:
    njs = Path(str(prefix) + ".njs")
    if not njs.is_file():
        return False, "BLAST .njs sidecar missing"
    try:
        meta = json.loads(njs.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return False, str(exc)
    nseq = meta.get("number-of-sequences")
    title = str(meta.get("description") or meta.get("title") or "")
    expected = GENCODE_V45_TRANSCRIPTS["n_transcripts"]
    if nseq != expected:
        return False, f"nseq={nseq} (expected {expected}) title={title}"
    return True, f"nseq={nseq} {title}".strip()
