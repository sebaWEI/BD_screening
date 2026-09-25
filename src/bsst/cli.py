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
    DBSNP_B151_GRCH38P7_COMMON_ALL,
    GENCODE_V45_TRANSCRIPTS,
    catalog,
    discover_blast_db,
    discover_variant_vcf,
)

app = typer.Typer(help="bsst: BLAST filter for Hepha antisense sites. RNAup ranking is off unless --rnaup.")
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
        None, help="Expected VCF SHA-256. Defaults to the pinned digest for --dbsnp-common-all."
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
            "Use NCBI dbSNP b151 GRCh38.p7 common_all_20180418 "
            f"({DBSNP_B151_GRCH38P7_COMMON_ALL['filename']})."
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
    console.print("BLAST+ and RNAup are external system tools and are not installed by uv.")


@app.command("check_requirements")
def check_requirements() -> None:
    """Report tools and configured data. Exits 1 if a required binary is missing or a pinned file does not match."""
    report = tool_report()
    table = Table("Resource", "Value", "Status")
    for row in report:
        table.add_row(row["resource"], row["value"], row["status"])
    console.print(table)
    console.print("A BLAST filter needs blastn, makeblastdb, blast_db, and gencode_utr_index `ok`.")
    console.print("RNAup is optional and used only with `bsst filter --rnaup`.")
    console.print("variant_vcf is optional and used only with `bsst filter --variants`.")
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
            + ". Finish `db init` before `bsst filter`, or pass `--variants` / `--rnaup` only when those tools are present.[/yellow]"
        )


@app.command("resources")
def resources_cmd() -> None:
    """Print the pinned scientific resources (producer names, not nicknames)."""
    table = Table("Role", "Producer / name", "Release", "Assembly", "URL / path")
    for item in catalog():
        table.add_row(
            str(item.get("role", "")),
            str(item.get("producer") or item.get("id")),
            str(
                item.get("variant_release")
                or item.get("transcriptome_release")
                or item.get("release")
                or item.get("ensembl_version")
                or item.get("task")
                or item.get("program")
                or ("live" if item.get("pinned") is False else "")
            ),
            str(item.get("assembly") or "—"),
            str(item.get("url") or item.get("base_url") or item.get("path") or item.get("manual") or item.get("docs_url") or ""),
        )
    console.print(table)
    console.print("Verification steps: docs/resources.md")


@config_app.command("show")
def config_show() -> None:
    console.print_json(json.dumps(load_config()))


def _options(
    window_size: int,
    step: int,
    context: int,
    temperature: float,
    include_both: bool,
    min_anchor_overlap: float,
    min_af: float | None,
    offtarget_min_length: int,
) -> SelectOptions:
    return SelectOptions(
        window_size=window_size,
        step=step,
        context=context,
        temperature=temperature,
        include_both=include_both,
        min_anchor_overlap=min_anchor_overlap,
        min_af=min_af,
        offtarget_min_length=offtarget_min_length,
    )


def _load_utr(
    fasta: Path,
    gene: str,
    chrom: str | None,
    start: int | None,
    end: int | None,
    strand: str | None,
) -> Target:
    if strand is not None and strand not in {"+", "-"}:
        raise typer.BadParameter("strand must be + or -")
    record = SeqIO.read(fasta, "fasta")
    header = _header_metadata(record.description)
    chrom = chrom or header.get("chrom")
    start = start if start is not None else (
        int(header["start"]) if "start" in header else None
    )
    end = end if end is not None else (
        int(header["end"]) if "end" in header else None
    )
    strand = strand or (
        header["strand"] if header.get("strand") in {"+", "-"} else "+"
    )
    transcript_id = None
    name_match = re.fullmatch(
        r"([A-Za-z0-9.-]+)_(ENST[0-9]+(?:\.[0-9]+)?)", record.id
    )
    if name_match:
        transcript_id = name_match.group(2)
    return Target(
        sequence=str(record.seq).upper(), name=record.id, chrom=chrom,
        start=start, end=end, strand=strand, gene=gene,
        transcript_id=transcript_id, source=f"FASTA: {fasta}",
    )


def _load_sites(fasta: Path) -> list[tuple[str, str]]:
    sites = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta, "fasta")]
    if not sites:
        raise typer.BadParameter("sites FASTA has no records")
    return sites


@app.command("filter")
def filter_sites(
    gene: str = typer.Option(..., "--gene", help="Gene symbol used to ignore self BLAST hits."),
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
    variant_vcf: Path | None = typer.Option(None, "--variant-vcf"),
    blast_db: str | None = typer.Option(None),
    variants: bool = typer.Option(False, "--variants", help="Drop sites overlapping dbSNP common_all."),
    rnaup: bool = typer.Option(False, "--rnaup", help="Score sites that passed the filter and rank by ΔG."),
    window_size: int = typer.Option(40),
    step: int = typer.Option(1),
    context: int = typer.Option(120),
    temperature: float = typer.Option(37.0),
    include_both: bool = typer.Option(False),
    min_anchor_overlap: float = typer.Option(1.0),
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
    resolved_vcf = variant_vcf or (
        Path(cfg["variant_vcf"]) if cfg.get("variant_vcf") else discover_variant_vcf()
    )
    resolved_blast = blast_db or cfg.get("blast_db") or discover_blast_db()
    run_dir = select(
        target,
        _options(
            window_size, step, context, temperature, include_both,
            min_anchor_overlap, min_af, offtarget_min_length,
        ),
        runs_dir=runs_dir or Path(cfg["runs_dir"]),
        variant_vcf=resolved_vcf,
        blast_db=resolved_blast,
        sites=_load_sites(sites) if sites else None,
        use_variants=variants,
        use_rnaup=rnaup,
        self_tokens=[token for token in (gene, target.transcript_id) if token],
    )
    console.print(f"Run completed: {run_dir}")


if __name__ == "__main__":
    app()
