from __future__ import annotations

import json
import re
from pathlib import Path

import typer
from Bio import SeqIO
from rich.console import Console
from rich.table import Table

from .config import load_config
from .db import initialize, tool_report
from .fetch import EnsemblArchiveError, fetch_gene_utr
from .models import SelectOptions, Target
from .pipeline import select
from .resources import (
    DBSNP_B157_GRCH38P14,
    GENCODE_V45_TRANSCRIPTS,
    catalog,
    discover_blast_db,
    discover_variant_vcf,
)

app = typer.Typer(
    help="Binding Site Safety Filter (bssf): BLAST + optional variant filter for Hepha antisense sites."
)
db_app = typer.Typer(help="Manage local database configuration.")
config_app = typer.Typer(help="Inspect configuration.")
app.add_typer(db_app, name="db")
app.add_typer(config_app, name="config")
console = Console()


def _header_metadata(description: str) -> dict[str, str]:
    """Parse the key=value header emitted by this project and older fetchers."""
    return dict(re.findall(r"\b(chrom|start|end|strand)=([^\s]+)", description))


@db_app.command("init")
def db_init(
    db_dir: Path | None = typer.Option(None, help="Database directory."),
    variant_vcf: Path | None = typer.Option(
        None, exists=True, readable=True, help="Configure an existing .vcf or .vcf.gz."
    ),
    blast_db: str | None = typer.Option(None, help="Configure an existing BLAST DB prefix."),
    variant_url: str | None = typer.Option(None, help="Download a VCF from this URL."),
    transcriptome_url: str | None = typer.Option(
        None, help="Download FASTA(.gz) and build a BLAST database."
    ),
    variant_sha256: str | None = typer.Option(
        None,
        help=(
            "Expected SHA-256 for a custom --variant-url download. "
            "Pinned --dbsnp-common-all verifies the NCBI source by MD5 instead."
        ),
    ),
    transcriptome_sha256: str | None = typer.Option(
        None,
        help="Expected downloaded FASTA(.gz) SHA-256. Defaults to the pinned archive digest for --gencode-v45-transcripts.",
    ),
    assembly: str | None = typer.Option(None, help="Reference assembly, e.g. GRCh38."),
    variant_source: str | None = typer.Option(None, help="Variant resource name."),
    variant_release: str | None = typer.Option(None, help="Variant resource release."),
    transcriptome_release: str | None = typer.Option(None, help="Transcript annotation release."),
    transcriptome_source: str | None = typer.Option(None, help="Transcriptome resource name."),
    transcriptome_assembly: str | None = typer.Option(
        None, help="Transcriptome genome assembly, e.g. GRCh38.p14."
    ),
    dbsnp_common_all: bool = typer.Option(
        False,
        "--dbsnp-common-all",
        help=(
            "Pin NCBI dbSNP b157 GRCh38.p14 for --variants: download the ~3 MB "
            f"tabix index ({DBSNP_B157_GRCH38P14['tbi_filename']}) and query "
            "INFO/COMMON sites on demand via HTTPS (no 28 GB VCF download)."
        ),
    ),
    gencode_v45_transcripts: bool = typer.Option(
        False,
        "--gencode-v45-transcripts",
        help=(
            "Use GENCODE 45 CHR transcripts "
            f"({GENCODE_V45_TRANSCRIPTS['filename']}) as the BLAST subject."
        ),
    ),
    retries: int = typer.Option(3, min=1, max=10),
    force: bool = typer.Option(False, help="Replace existing downloaded data."),
) -> None:
    directory, path = initialize(
        db_dir,
        variant_vcf=variant_vcf,
        blast_db=blast_db,
        variant_url=variant_url,
        transcriptome_url=transcriptome_url,
        variant_sha256=variant_sha256,
        transcriptome_sha256=transcriptome_sha256,
        assembly=assembly,
        variant_source=variant_source,
        variant_release=variant_release,
        transcriptome_release=transcriptome_release,
        transcriptome_source=transcriptome_source,
        transcriptome_assembly=transcriptome_assembly,
        dbsnp_common_all=dbsnp_common_all,
        gencode_v45_transcripts=gencode_v45_transcripts,
        retries=retries,
        force=force,
    )
    console.print(f"Database directory: {directory}")
    console.print(f"Configuration: {path}")
    console.print("BLAST+ is an external system tool and is not installed by uv.")


@app.command("check_requirements")
def check_requirements() -> None:
    """Report tools and configured data. Exits 1 if a required binary is missing or a pinned file does not match."""
    report = tool_report()
    table = Table("Resource", "Value", "Status")
    for row in report:
        table.add_row(row["resource"], row["value"], row["status"])
    console.print(table)
    console.print("A BLAST filter needs blastn, makeblastdb, blast_db, and gencode_utr_index `ok`.")
    console.print("variant_vcf is optional and used only with `bssf filter --variants`.")
    console.print("tabix is optional but recommended for large VCF region queries.")
    failed = [row for row in report if row["status"] in {"missing", "mismatch"}]
    if failed:
        names = ", ".join(row["resource"] for row in failed)
        console.print(f"[red]Requirements incomplete: {names}[/red]")
        raise typer.Exit(code=1)
    optional_missing = [row["resource"] for row in report if row["status"] == "missing/optional"]
    if optional_missing:
        console.print(
            "[yellow]Optional items still missing: "
            + ", ".join(optional_missing)
            + ". Finish `db init` before `bssf filter`, or pass `--variants` only when the VCF is present.[/yellow]"
        )


@app.command("resources")
def resources_cmd() -> None:
    """Print the pinned scientific resources (producer names, not nicknames)."""
    table = Table("Role", "Producer / name", "Release", "Assembly", "URL / path")
    for item in catalog():
        table.add_row(
            str(item.get("role", "")),
            str(item.get("producer") or item.get("program") or item.get("id", "")),
            str(item.get("release") or item.get("ensembl_version") or ""),
            str(item.get("assembly") or ""),
            str(item.get("url") or item.get("path") or item.get("manual") or ""),
        )
    console.print(table)


@config_app.command("show")
def config_show() -> None:
    console.print_json(json.dumps(load_config()))


def _options(
    window_size: int,
    step: int,
    min_af: float | None,
    offtarget_min_length: int,
) -> SelectOptions:
    return SelectOptions(
        window_size=window_size,
        step=step,
        min_af=min_af,
        offtarget_min_length=offtarget_min_length,
    )


def _load_utr(
    utr: Path,
    gene: str,
    chrom: str | None,
    start: int | None,
    end: int | None,
    strand: str | None,
) -> Target:
    records = list(SeqIO.parse(str(utr), "fasta"))
    if len(records) != 1:
        raise typer.BadParameter(f"--utr must contain exactly one sequence, found {len(records)}")
    record = records[0]
    meta = _header_metadata(record.description)
    resolved_strand = strand or meta.get("strand") or "+"
    resolved_chrom = chrom or meta.get("chrom")
    resolved_start = start if start is not None else (
        int(meta["start"]) if "start" in meta else None
    )
    resolved_end = end if end is not None else (
        int(meta["end"]) if "end" in meta else None
    )
    return Target(
        str(record.seq),
        record.id,
        chrom=resolved_chrom,
        start=resolved_start,
        end=resolved_end,
        strand=resolved_strand,
        gene=gene,
    )


def _load_sites(path: Path) -> list[tuple[str, str]]:
    return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(str(path), "fasta")]


@app.command("filter")
def filter_cmd(
    gene: str = typer.Option(..., "--gene", help="Gene symbol; required so self BLAST hits are ignored."),
    utr: Path | None = typer.Option(
        None, "--utr", exists=True, readable=True,
        help="3'UTR FASTA. Omitted sequences are fetched from Ensembl 111.",
    ),
    sites: Path | None = typer.Option(
        None, "--sites", exists=True, readable=True,
        help="Sense binding sites. Omitted sites are 40 nt windows stepped by 1.",
    ),
    runs_dir: Path | None = typer.Option(None, help="Override configured runs directory."),
    species: str = typer.Option("homo_sapiens"),
    chrom: str | None = typer.Option(None),
    start: int | None = typer.Option(None, help="0-based BED start of the UTR."),
    end: int | None = typer.Option(None, help="0-based BED end of the UTR."),
    strand: str | None = typer.Option(None, help="+ or -; inferred from a compatible UTR header, otherwise +."),
    variant_vcf: str | None = typer.Option(
        None, "--variant-vcf", help="Local VCF(.gz) path or HTTPS bgzip VCF URL for tabix."
    ),
    blast_db: str | None = typer.Option(None),
    variants: bool = typer.Option(
        False, "--variants", help="Drop sites overlapping dbSNP COMMON (tabix region query)."
    ),
    window_size: int = typer.Option(40),
    step: int = typer.Option(1),
    min_af: float | None = typer.Option(
        None, min=0.0, max=1.0, help="With --variants, ignore variants below this frequency."
    ),
    offtarget_min_length: int = typer.Option(
        20, "--offtarget-min-length", min=1,
        help="Fail a site when a non-self BLAST alignment is at least this many nt. Default 20.",
    ),
) -> None:
    if utr is None:
        try:
            target = fetch_gene_utr(gene, species)
        except EnsemblArchiveError as exc:
            console.print(f"[red]{exc}[/red]")
            raise typer.Exit(code=1) from exc
    else:
        target = _load_utr(utr, gene, chrom, start, end, strand)
    if variants and (target.chrom is None or target.start is None or target.end is None):
        console.print("[yellow]Warning: UTR has no complete coordinates; variant filtering will be skipped.[/yellow]")
    cfg = load_config()
    resolved_vcf = variant_vcf or cfg.get("variant_vcf") or discover_variant_vcf()
    resolved_blast = blast_db or cfg.get("blast_db") or discover_blast_db()
    run_dir = select(
        target,
        _options(window_size, step, min_af, offtarget_min_length),
        runs_dir=runs_dir or Path(cfg["runs_dir"]),
        variant_vcf=resolved_vcf,
        blast_db=resolved_blast,
        sites=_load_sites(sites) if sites else None,
        use_variants=variants,
        self_tokens=[token for token in (gene, target.transcript_id) if token],
    )
    console.print(f"Run completed: {run_dir}")


if __name__ == "__main__":
    app()
