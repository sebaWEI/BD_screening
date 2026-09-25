from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path
from typing import Any

from .config import load_config, save_config
from .pipeline import resolve_executable
from .resources import (
    DBSNP_B151_GRCH38P7_COMMON_ALL,
    GENCODE_V45_TRANSCRIPTS,
    blast_db_is_present,
    bundled_blast_db_prefix,
    bundled_dbsnp_vcf,
    bundled_gencode_fasta,
    data_dir,
    discover_blast_db,
    discover_variant_vcf,
    tabix_index_path,
    verify_blast_db,
    verify_dbsnp_vcf,
    verify_gencode_fasta,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(
    url: str,
    destination: Path,
    *,
    expected_sha256: str | None = None,
    retries: int = 3,
    force: bool = False,
) -> Path:
    """Stream a URL to disk with bounded retries and optional SHA-256 verification."""
    if destination.exists() and not force:
        if expected_sha256 and _sha256(destination) != expected_sha256.lower():
            raise ValueError(f"checksum mismatch for existing file: {destination}")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "bsst/0.3"})
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if expected_sha256 and _sha256(partial) != expected_sha256.lower():
                raise ValueError(f"checksum mismatch for downloaded file: {url}")
            partial.replace(destination)
            return destination
        except ValueError:
            partial.unlink(missing_ok=True)
            raise
        except Exception as exc:
            last_error = exc
            partial.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 8))
    raise RuntimeError(f"download failed after {retries} attempts: {url}") from last_error


def initialize(
    db_dir: Path | None = None,
    *,
    variant_vcf: Path | None = None,
    blast_db: str | None = None,
    variant_url: str | None = None,
    transcriptome_url: str | None = None,
    variant_sha256: str | None = None,
    transcriptome_sha256: str | None = None,
    assembly: str | None = None,
    variant_source: str | None = None,
    variant_release: str | None = None,
    transcriptome_release: str | None = None,
    transcriptome_source: str | None = None,
    transcriptome_assembly: str | None = None,
    dbsnp_common_all: bool = False,
    gencode_v45_transcripts: bool = False,
    retries: int = 3,
    force: bool = False,
) -> tuple[Path, Path]:
    """Create/configure data resources and optionally download/build them."""
    cfg = load_config()
    directory = db_dir or Path(cfg["db_dir"])
    directory.mkdir(parents=True, exist_ok=True)
    cfg["db_dir"] = str(directory.resolve())
    spec = DBSNP_B151_GRCH38P7_COMMON_ALL
    gencode = GENCODE_V45_TRANSCRIPTS

    if dbsnp_common_all:
        assembly = assembly or spec["assembly"]
        variant_source = variant_source or spec["variant_source"]
        variant_release = variant_release or spec["variant_release"]
        destination = bundled_dbsnp_vcf()
        if variant_vcf is None:
            variant_vcf = download_file(
                variant_url or spec["url"],
                destination,
                expected_sha256=variant_sha256 or spec["sha256"],
                retries=retries,
                force=force,
            )
        variant_url = None

    if variant_url:
        suffix = ".vcf.gz" if variant_url.lower().endswith(".gz") else ".vcf"
        variant_vcf = download_file(
            variant_url,
            directory / f"variants{suffix}",
            expected_sha256=variant_sha256,
            retries=retries,
            force=force,
        )
    if variant_vcf:
        if not variant_vcf.exists():
            raise FileNotFoundError(variant_vcf)
        cfg["variant_vcf"] = str(variant_vcf.resolve())
        ensure_tabix(variant_vcf)

    if gencode_v45_transcripts:
        transcriptome_source = transcriptome_source or gencode["transcriptome_source"]
        transcriptome_release = transcriptome_release or gencode["transcriptome_release"]
        transcriptome_assembly = transcriptome_assembly or gencode["assembly"]
        fasta = bundled_gencode_fasta()
        archive_digest = transcriptome_sha256 or gencode["archive_sha256"]
        fasta_digest = None if transcriptome_sha256 else gencode["fasta_sha256"]
        if transcriptome_url is None and (not fasta.is_file() or force):
            archive = download_file(
                gencode["url"],
                data_dir() / gencode["archive_filename"],
                expected_sha256=archive_digest,
                retries=retries,
                force=force,
            )
            if force or not fasta.is_file():
                fasta.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(archive, "rb") as source, fasta.open("wb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
        if transcriptome_url is None and fasta.is_file() and fasta_digest:
            if _sha256(fasta) != fasta_digest.lower():
                raise ValueError(f"checksum mismatch for existing file: {fasta}")
        if blast_db is None:
            prefix = bundled_blast_db_prefix()
            if force or not blast_db_is_present(prefix):
                _make_blast_db(fasta, prefix, title=gencode["blast_title"])
            blast_db = str(prefix)
        transcriptome_url = None

    if transcriptome_url:
        archive = download_file(
            transcriptome_url,
            directory / (
                "transcriptome.fa.gz"
                if transcriptome_url.lower().endswith(".gz")
                else "transcriptome.fa"
            ),
            expected_sha256=transcriptome_sha256,
            retries=retries,
            force=force,
        )
        fasta = directory / "transcriptome.fa"
        if archive.name.endswith(".gz"):
            if force or not fasta.exists():
                with gzip.open(archive, "rb") as source, fasta.open("wb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
        else:
            fasta = archive
        prefix = directory / "transcriptome_db"
        _make_blast_db(fasta, prefix)
        blast_db = str(prefix)
    if blast_db:
        cfg["blast_db"] = str(Path(blast_db).expanduser())
    for key, value in {
        "assembly": assembly,
        "variant_source": variant_source,
        "variant_release": variant_release,
        "transcriptome_release": transcriptome_release,
        "transcriptome_source": transcriptome_source,
        "transcriptome_assembly": transcriptome_assembly,
    }.items():
        if value:
            cfg[key] = value

    return directory, save_config(cfg)


def ensure_tabix(vcf: Path) -> Path | None:
    """Build a tabix index for a bgzip VCF when the tabix binary is present."""
    if not vcf.name.endswith(".vcf.gz"):
        return None
    index = tabix_index_path(vcf)
    if index.is_file():
        return index
    tabix = shutil.which("tabix")
    if not tabix:
        return None
    subprocess.run([tabix, "-p", "vcf", str(vcf)], check=True)
    return index if index.is_file() else None


def _row(
    name: str,
    value: Any,
    ok: bool,
    *,
    missing_ok: bool = True,
    required: bool = False,
) -> dict[str, str]:
    if ok:
        status = "ok"
    elif required:
        status = "missing"
    elif missing_ok:
        status = "missing/optional"
    else:
        status = "mismatch"
    return {"resource": name, "value": str(value or "not found"), "status": status}


def tool_report() -> list[dict[str, str]]:
    cfg = load_config()
    blast_db = cfg.get("blast_db") or discover_blast_db()
    variant_vcf = cfg.get("variant_vcf") or (
        str(discover_variant_vcf()) if discover_variant_vcf() else None
    )
    vcf_path = Path(variant_vcf) if variant_vcf else None
    fasta = bundled_gencode_fasta()
    rows = [
        _row(
            "RNAup",
            resolve_executable("RNAup"),
            bool(resolve_executable("RNAup")),
        ),
        _row(
            "blastn",
            resolve_executable("blastn"),
            bool(resolve_executable("blastn")),
            required=True,
        ),
        _row(
            "makeblastdb",
            shutil.which("makeblastdb"),
            bool(shutil.which("makeblastdb")),
            required=True,
        ),
        _row("tabix", shutil.which("tabix"), bool(shutil.which("tabix"))),
    ]
    blast_ok = bool(blast_db and blast_db_is_present(blast_db))
    rows.append(_row("blast_db", blast_db, blast_ok))
    if blast_ok:
        identity_ok, detail = verify_blast_db(blast_db)
        rows.append(_row("blast_db_identity", detail, identity_ok, missing_ok=False))
    if fasta.is_file():
        fasta_ok, fasta_detail = verify_gencode_fasta(fasta)
        rows.append(_row("gencode_fasta", fasta, True))
        rows.append(_row("gencode_fasta_identity", fasta_detail, fasta_ok, missing_ok=False))
    else:
        rows.append(_row("gencode_fasta", None, False))
    rows.extend(
        [
            _row("transcriptome_source", cfg.get("transcriptome_source"), bool(cfg.get("transcriptome_source"))),
            _row("transcriptome_release", cfg.get("transcriptome_release"), bool(cfg.get("transcriptome_release"))),
            _row("transcriptome_assembly", cfg.get("transcriptome_assembly"), bool(cfg.get("transcriptome_assembly"))),
        ]
    )
    vcf_exists = bool(vcf_path and vcf_path.exists())
    rows.append(_row("variant_vcf", variant_vcf, vcf_exists))
    if vcf_exists and vcf_path is not None:
        header_ok, header_detail = verify_dbsnp_vcf(vcf_path)
        rows.append(_row("variant_vcf_identity", header_detail, header_ok, missing_ok=False))
        tbi = tabix_index_path(vcf_path)
        rows.append(_row("variant_tabix", tbi if tbi.is_file() else None, tbi.is_file()))
    rows.extend(
        [
            _row("variant_assembly", cfg.get("assembly"), bool(cfg.get("assembly"))),
            _row("variant_source", cfg.get("variant_source"), bool(cfg.get("variant_source"))),
            _row("variant_release", cfg.get("variant_release"), bool(cfg.get("variant_release"))),
        ]
    )
    return rows


def _make_blast_db(fasta: Path, prefix: Path, *, title: str | None = None) -> None:
    if not fasta.is_file():
        raise FileNotFoundError(fasta)
    makeblastdb = shutil.which("makeblastdb")
    if not makeblastdb:
        raise RuntimeError("makeblastdb is required to build the transcriptome database")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [makeblastdb, "-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix)]
    if title:
        command.extend(["-title", title])
    subprocess.run(command, check=True)
