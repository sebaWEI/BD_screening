"""Pinned scientific resources used by bssf.

Names, URLs, and version strings are the identifiers assigned by the data
producer. Do not rename a file's contents after another database.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# NCBI dbSNP (GRCh38.p14) — tabix region queries
#
# NCBI no longer ships a separate common_all VCF on p14. The pinned source is
# GCF_000001405.40 (~28 GB). bssf does **not** download that file: ``db init
# --dbsnp-common-all`` fetches only the ~3 MB ``.tbi``, and ``--variants``
# runs ``tabix`` over the HTTPS URL for the target UTR interval, keeping
# INFO/COMMON records and remapping NC_* → 1..22/X/Y/MT.
# ---------------------------------------------------------------------------
DBSNP_B157_GRCH38P14: dict[str, Any] = {
    "id": "dbsnp_b157_grch38p14_tabix",
    "role": "variant_vcf",
    "producer": "NCBI dbSNP",
    "source_filename": "GCF_000001405.40.gz",
    "tbi_filename": "GCF_000001405.40.gz.tbi",
    # Prefer ftp.ncbi.nlm.nih.gov (often faster / more reachable than ftp.ncbi.nih.gov).
    "url": "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
    "directory_url": "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/",
    "tbi_url": "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi",
    "url_alt": "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
    "assembly": "GRCh38.p14",
    "refseq_accession": "GCF_000001405.40",
    "variant_source": "NCBI dbSNP",
    "variant_release": "b157 fileDate=20241205 tabix INFO/COMMON",
    "chrom_style": "unprefixed",  # after remap: 1, 2, ..., X, Y, MT — not chr1
    "source_chrom_style": "refseq",  # NC_000001.11 …
    "frequency_fields": ("AF", "CAF", "FREQ"),
    "common_flag": "COMMON",
    "common_only_default": True,
    "header_checks": (
        "##reference=GRCh38.p14",
        "##dbSNP_BUILD_ID=157",
    ),
    # Primary assembled molecules only (GRCh38.p14 assembly report).
    "refseq_to_chrom": {
        "NC_000001.11": "1",
        "NC_000002.12": "2",
        "NC_000003.12": "3",
        "NC_000004.12": "4",
        "NC_000005.10": "5",
        "NC_000006.12": "6",
        "NC_000007.14": "7",
        "NC_000008.11": "8",
        "NC_000009.12": "9",
        "NC_000010.11": "10",
        "NC_000011.10": "11",
        "NC_000012.12": "12",
        "NC_000013.11": "13",
        "NC_000014.9": "14",
        "NC_000015.10": "15",
        "NC_000016.10": "16",
        "NC_000017.11": "17",
        "NC_000018.10": "18",
        "NC_000019.10": "19",
        "NC_000020.11": "20",
        "NC_000021.9": "21",
        "NC_000022.11": "22",
        "NC_000023.11": "X",
        "NC_000024.10": "Y",
        "NC_012920.1": "MT",
    },
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
    "producer": "bssf examples (Ensembl-style 3′UTR FASTA)",
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
    "producer": "bssf examples (Ensembl-style 3′UTR FASTA)",
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
        DBSNP_B157_GRCH38P14,
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
    if (repo / "pyproject.toml").is_file() and (repo / "src" / "bssf").is_dir():
        return repo
    return Path.cwd()


def data_dir() -> Path:
    return package_root() / "data"


def bundled_dbsnp_source() -> Path:
    return data_dir() / DBSNP_B157_GRCH38P14["source_filename"]


def bundled_dbsnp_tbi() -> Path:
    return data_dir() / DBSNP_B157_GRCH38P14["tbi_filename"]


def chrom_to_refseq(chrom: str) -> str | None:
    """Map unprefixed / chr-prefixed names to GCF RefSeq contig IDs."""
    key = str(chrom).strip()
    if key.lower().startswith("chr"):
        key = key[3:]
    upper = key.upper()
    if upper == "M":
        upper = "MT"
    for contig, name in DBSNP_B157_GRCH38P14["refseq_to_chrom"].items():
        if name.upper() == upper:
            return contig
    return None


def refseq_to_chrom(contig: str) -> str:
    return DBSNP_B157_GRCH38P14["refseq_to_chrom"].get(contig, contig)


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


def discover_variant_vcf() -> str | None:
    """Pinned HTTPS GCF URL once the local ``.tbi`` has been fetched."""
    if bundled_dbsnp_tbi().is_file():
        return DBSNP_B157_GRCH38P14["url"]
    return None


def discover_blast_db() -> str | None:
    prefix = bundled_blast_db_prefix()
    return str(prefix) if blast_db_is_present(prefix) else None


def tabix_index_path(vcf: Path | str) -> Path:
    return Path(str(vcf) + ".tbi")


def verify_dbsnp_tabix(source: str | Path) -> tuple[bool, str]:
    """Confirm tabix can read the pinned GRCh38.p14 / b157 header."""
    import shutil
    import subprocess

    tabix = shutil.which("tabix")
    if not tabix:
        return False, "tabix not on PATH"
    try:
        proc = subprocess.run(
            [tabix, "-H", str(source)],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return False, str(exc)
    if proc.returncode != 0:
        return False, (proc.stderr or proc.stdout or "tabix -H failed").strip()[:200]
    text = proc.stdout
    missing = [
        field for field in DBSNP_B157_GRCH38P14["header_checks"]
        if field not in text
    ]
    if missing:
        return False, "missing " + ", ".join(missing)
    return True, "dbSNP b157 GRCh38.p14 tabix"


# Back-compat alias used by older call sites.
def verify_dbsnp_vcf(path: Path | str) -> tuple[bool, str]:
    return verify_dbsnp_tabix(path)


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
