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
from .annotation import write_utr_index
from .resources import (
    DBSNP_B157_GRCH38P14,
    GENCODE_V45_ANNOTATION,
    GENCODE_V45_TRANSCRIPTS,
    blast_db_is_present,
    bundled_blast_db_prefix,
    bundled_dbsnp_tbi,
    bundled_gencode_fasta,
    bundled_gencode_gff3,
    bundled_utr_index,
    data_dir,
    discover_blast_db,
    discover_variant_vcf,
    tabix_index_path,
    verify_blast_db,
    verify_dbsnp_tabix,
    verify_gencode_fasta,
)


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    return _digest(path, "sha256")


def download_file(
    url: str,
    destination: Path,
    *,
    expected_sha256: str | None = None,
    expected_md5: str | None = None,
    retries: int = 3,
    force: bool = False,
) -> Path:
    """Stream a URL to disk with bounded retries and optional digest verification."""
    if destination.exists() and not force:
        if expected_sha256 and _sha256(destination) != expected_sha256.lower():
            raise ValueError(f"checksum mismatch for existing file: {destination}")
        if expected_md5 and _digest(destination, "md5") != expected_md5.lower():
            raise ValueError(f"checksum mismatch for existing file: {destination}")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "bssf/0.4"})
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            if expected_sha256 and _sha256(partial) != expected_sha256.lower():
                raise ValueError(f"checksum mismatch for downloaded file: {url}")
            if expected_md5 and _digest(partial, "md5") != expected_md5.lower():
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
    spec = DBSNP_B157_GRCH38P14
    gencode = GENCODE_V45_TRANSCRIPTS

    if dbsnp_common_all:
        assembly = assembly or spec["assembly"]
        variant_source = variant_source or spec["variant_source"]
        variant_release = variant_release or spec["variant_release"]
        # Only the tabix index (~3 MB). Queries hit the HTTPS GCF URL on demand.
        download_file(
            spec["tbi_url"],
            bundled_dbsnp_tbi(),
            retries=retries,
            force=force,
        )
        cfg["variant_vcf"] = variant_url or spec["url"]
        variant_url = None
        variant_vcf = None

    if variant_url:
        if variant_url.startswith(("http://", "https://")):
            cfg["variant_vcf"] = variant_url
            variant_vcf = None
        else:
            suffix = ".vcf.gz" if variant_url.lower().endswith(".gz") else ".vcf"
            variant_vcf = download_file(
                variant_url,
                directory / f"variants{suffix}",
                expected_sha256=variant_sha256,
                retries=retries,
                force=force,
            )
    if variant_vcf:
        if not Path(variant_vcf).exists():
            raise FileNotFoundError(variant_vcf)
        cfg["variant_vcf"] = str(Path(variant_vcf).resolve())
        ensure_tabix(Path(variant_vcf))

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
        gff3 = bundled_gencode_gff3()
        expected_md5 = GENCODE_V45_ANNOTATION["archive_md5"]
        if gff3.is_file() and _digest(gff3, "md5") != expected_md5:
            gff3.unlink()
        if not gff3.is_file() or force:
            download_file(
                GENCODE_V45_ANNOTATION["url"],
                gff3,
                expected_sha256=GENCODE_V45_ANNOTATION["archive_sha256"] or None,
                retries=retries,
                force=force,
            )
        if _digest(gff3, "md5") != expected_md5:
            raise ValueError(f"checksum mismatch for existing file: {gff3}")
        index = bundled_utr_index()
        if force or not index.is_file():
            write_utr_index(gff3, index)
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
    variant_vcf = cfg.get("variant_vcf") or discover_variant_vcf()
    fasta = bundled_gencode_fasta()
    rows = [
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
        _row("tabix", shutil.which("tabix"), bool(shutil.which("tabix")), required=False),
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
    utr_index = bundled_utr_index()
    rows.append(_row("gencode_utr_index", utr_index if utr_index.is_file() else None, utr_index.is_file(), required=True))
    rows.extend(
        [
            _row("transcriptome_source", cfg.get("transcriptome_source"), bool(cfg.get("transcriptome_source"))),
            _row("transcriptome_release", cfg.get("transcriptome_release"), bool(cfg.get("transcriptome_release"))),
            _row("transcriptome_assembly", cfg.get("transcriptome_assembly"), bool(cfg.get("transcriptome_assembly"))),
        ]
    )
    http = bool(variant_vcf and str(variant_vcf).startswith(("http://", "https://")))
    local = Path(variant_vcf) if variant_vcf and not http else None
    vcf_ok = http or bool(local and local.exists())
    rows.append(_row("variant_vcf", variant_vcf, vcf_ok))
    if vcf_ok and variant_vcf:
        if shutil.which("tabix"):
            header_ok, header_detail = verify_dbsnp_tabix(variant_vcf)
            rows.append(
                _row("variant_vcf_identity", header_detail, header_ok, missing_ok=False)
            )
        else:
            rows.append(
                _row("variant_vcf_identity", "tabix not on PATH", False, missing_ok=True)
            )
        tbi = bundled_dbsnp_tbi() if http else (tabix_index_path(local) if local else None)
        rows.append(
            _row("variant_tabix", tbi if tbi and tbi.is_file() else None, bool(tbi and tbi.is_file()))
        )
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
