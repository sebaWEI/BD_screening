from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from uuid import uuid4

import pandas as pd
from Bio.Seq import Seq

from .config import load_config
from .logging import log_command, run_logger, write_manifest
from .models import MODEL_VERSION, SelectOptions, Target
from .annotation import load_utr_index, utr_regions
from .resources import BLASTN_SHORT, discover_utr_index

WINDOW_COLUMNS = [
    "name", "chrom", "start", "end", "strand", "utr_start", "utr_end",
    "target_sequence", "bd_sequence", "gc_fraction", "status", "failure_reason",
    "model_version",
]
BLAST_COLUMNS = [
    "qseqid", "sseqid", "stitle", "sallacc", "pident", "length", "qlen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def resolve_executable(name: str, explicit: str | None = None) -> str | None:
    if explicit:
        return explicit
    env_value = os.environ.get(f"BSSF_{name.upper()}")
    configured = load_config().get(f"{name.lower()}_exe")
    candidates = [
        env_value,
        configured,
        shutil.which(name),
        str(Path.cwd() / ".venv" / "bin" / name),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate))
    return None


def generate_windows(target: Target, window_size: int = 40, step: int = 1) -> pd.DataFrame:
    if window_size <= 0 or step <= 0:
        raise ValueError("window_size and step must be positive")
    seq = target.sequence.upper().replace("U", "T")
    rows = []
    for offset in range(0, len(seq) - window_size + 1, step):
        window = seq[offset : offset + window_size]
        rows.append(_window_row(target, f"{target.name}_w{offset}", window, offset))
    return _windows_frame(rows)


def _genomic_interval(target: Target, offset: int, length: int) -> tuple[object, object]:
    if target.start is None or target.end is None:
        return pd.NA, pd.NA
    if target.strand == "+":
        return target.start + offset, target.start + offset + length
    return target.end - offset - length, target.end - offset


def _window_row(
    target: Target,
    name: str,
    window: str,
    offset: int,
    *,
    status: str = "pending",
    failure_reason: str = "",
) -> dict[str, object]:
    start, end = _genomic_interval(target, offset, len(window))
    return {
        "name": name,
        "chrom": target.chrom,
        "start": start,
        "end": end,
        "strand": target.strand,
        "utr_start": offset,
        "utr_end": offset + len(window),
        "target_sequence": window,
        "bd_sequence": str(Seq(window).reverse_complement()),
        "gc_fraction": (window.count("G") + window.count("C")) / len(window) if window else pd.NA,
        "status": status,
        "failure_reason": failure_reason,
        "model_version": MODEL_VERSION,
    }


def _windows_frame(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=[
            "name", "chrom", "start", "end", "strand", "utr_start", "utr_end",
            "target_sequence", "bd_sequence", "gc_fraction", "status",
            "failure_reason", "model_version",
        ],
    )


def sites_on_utr(target: Target, sites: list[tuple[str, str]]) -> pd.DataFrame:
    """Place sense sites on the UTR. Zero or multiple matches are not BLASTed."""
    seq = target.sequence.upper().replace("U", "T")
    rows: list[dict[str, object]] = []
    for name, raw in sites:
        site = raw.upper().replace("U", "T")
        positions: list[int] = []
        cursor = 0
        while site and (found := seq.find(site, cursor)) >= 0:
            positions.append(found)
            cursor = found + 1
        if len(positions) != 1:
            reason = "ambiguous" if len(positions) > 1 else "not_in_utr"
            rows.append(
                _window_row(target, name, site, 0, status="failed", failure_reason=reason)
            )
            rows[-1]["utr_start"] = pd.NA
            rows[-1]["utr_end"] = pd.NA
            rows[-1]["start"] = pd.NA
            rows[-1]["end"] = pd.NA
            continue
        rows.append(_window_row(target, name, site, positions[0]))
    return _windows_frame(rows)


def complexity_failure(
    sequence: str,
    min_gc: float = 0.3,
    max_gc: float = 0.7,
    max_homopolymer: int = 4,
) -> str | None:
    seq = sequence.upper()
    if not seq:
        return "empty_sequence"
    if re.search(rf"(.)\1{{{max_homopolymer},}}", seq):
        return "homopolymer"
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    return None if min_gc <= gc <= max_gc else "gc_out_of_range"


def _chrom(value: str) -> str:
    normalized = str(value).strip().lower()
    return normalized[3:] if normalized.startswith("chr") else normalized


def _freq_field_max_alt(freq: str) -> float | None:
    """Max alternate frequency across dbSNP ``FREQ=study:ref,alt,…|…`` blocks."""
    best: float | None = None
    for block in freq.split("|"):
        if ":" not in block:
            continue
        _, values = block.split(":", 1)
        nums = [float(x) for x in values.split(",") if x not in {"", "."}]
        if len(nums) < 2:
            continue
        alt = max(nums[1:])
        best = alt if best is None else max(best, alt)
    return best


def _info_allele_frequency(info: dict[str, str]) -> float | None:
    """Return a scalar frequency for filtering, or None if none is declared.

    Preference order:
    1. VCF ``AF`` (alternate-allele frequencies)
    2. dbSNP ``CAF`` (first value reference, then alternates)
    3. dbSNP ``FREQ`` (per-study ref,alt,… lists)
    """
    af_values = [
        float(x) for x in info.get("AF", "").split(",")
        if x not in {"", "."}
    ]
    if af_values:
        return max(af_values)
    caf_values = [
        float(x) for x in info.get("CAF", "").split(",")
        if x not in {"", "."}
    ]
    if len(caf_values) >= 2:
        return max(caf_values[1:])
    if info.get("FREQ"):
        return _freq_field_max_alt(info["FREQ"])
    return None


def _info_has_common(info: str) -> bool:
    return any(part == "COMMON" or part.startswith("COMMON=") for part in info.split(";"))


def _parse_vcf_record(
    line: str,
    min_af: float | None,
    *,
    common_only: bool = False,
) -> tuple[str, int, int, float | None] | None:
    from .resources import refseq_to_chrom

    if not line.strip() or line.startswith("#"):
        return None
    fields = line.rstrip().split("\t")
    if len(fields) < 8:
        return None
    if common_only and not _info_has_common(fields[7]):
        return None
    pos0 = int(fields[1]) - 1
    ref = fields[3]
    info = dict(
        item.split("=", 1) if "=" in item else (item, "")
        for item in fields[7].split(";")
    )
    af = _info_allele_frequency(info)
    if min_af is not None and (af is None or af < min_af):
        return None
    chrom = refseq_to_chrom(fields[0])
    return (_chrom(chrom), pos0, pos0 + max(1, len(ref)), af)


def _record_in_region(
    record: tuple[str, int, int, float | None],
    chrom: str | None,
    start: int | None,
    end: int | None,
) -> bool:
    if chrom is None:
        return True
    if record[0] != _chrom(chrom):
        return False
    if start is None or end is None:
        return True
    return start < record[2] and record[1] < end


def _is_http_source(source: str | Path) -> bool:
    return str(source).startswith(("http://", "https://"))


def _tabix_query_contig(chrom: str, *, gcf: bool) -> str:
    """Contig ID to pass to tabix (RefSeq for GCF, else normalized chrom)."""
    if gcf:
        from .resources import chrom_to_refseq

        return chrom_to_refseq(chrom) or _chrom(chrom)
    return _chrom(chrom)


def _tabix_vcf_lines(
    source: str | Path, chrom: str, start: int, end: int, *, gcf: bool = False
) -> list[str] | None:
    tabix = shutil.which("tabix")
    if not tabix:
        return None
    source_s = str(source)
    if not _is_http_source(source_s):
        index = Path(source_s + ".tbi")
        if not index.is_file():
            return None
    contig = _tabix_query_contig(chrom, gcf=gcf)
    region = f"{contig}:{start + 1}-{end}"
    proc = subprocess.run(
        [tabix, source_s, region],
        capture_output=True,
        text=True,
        timeout=180,
    )
    if proc.returncode != 0:
        return None
    return [line for line in proc.stdout.splitlines() if line.strip()]


def _looks_like_gcf_dbsnp(source: str | Path) -> bool:
    name = str(source)
    return "GCF_000001405.40" in name or name.endswith("GCF_000001405.40.gz")


def read_variant_vcf(
    path: str | Path,
    min_af: float | None = None,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    common_only: bool | None = None,
) -> list[tuple[str, int, int, float | None]]:
    """Read variants from a local VCF(.gz) or an HTTPS bgzip VCF via tabix.

    For the pinned NCBI GCF dbSNP release, ``common_only`` defaults to True and
    a genomic region is required (full-file scans of ~28 GB are not supported).
    """
    source = str(path)
    gcf = _looks_like_gcf_dbsnp(source)
    if common_only is None:
        common_only = gcf
    if gcf and (chrom is None or start is None or end is None):
        raise ValueError("GCF dbSNP queries require chrom/start/end (tabix region)")
    if not _is_http_source(source):
        path_obj = Path(source)
        if path_obj.suffix not in {".vcf", ".gz"} or (
            path_obj.suffix == ".gz" and not path_obj.name.endswith((".vcf.gz", ".40.gz"))
        ):
            # Allow GCF_000001405.40.gz (no .vcf in the name).
            if not path_obj.name.endswith(".gz") and path_obj.suffix != ".vcf":
                raise ValueError("variant file must end in .vcf or .vcf.gz")
    region = chrom is not None and start is not None and end is not None
    if region:
        tabix_lines = _tabix_vcf_lines(source, chrom, start, end, gcf=gcf)
        if tabix_lines is not None:
            records = []
            for line in tabix_lines:
                parsed = _parse_vcf_record(line, min_af, common_only=common_only)
                if parsed and _record_in_region(parsed, chrom, start, end):
                    records.append(parsed)
            return records
        if _is_http_source(source) or gcf:
            raise RuntimeError(
                f"tabix region query failed for {source} "
                f"({_tabix_query_contig(chrom, gcf=gcf)}:{start + 1}-{end}); "
                "is tabix installed and the network reachable?"
            )
    if _is_http_source(source):
        raise ValueError("HTTP variant sources require chrom/start/end for tabix")
    path_obj = Path(source)
    opener = gzip.open if path_obj.name.endswith(".gz") else open
    variants = []
    with opener(path_obj, "rt", encoding="utf-8") as handle:
        for line in handle:
            parsed = _parse_vcf_record(line, min_af, common_only=common_only)
            if parsed and _record_in_region(parsed, chrom, start, end):
                variants.append(parsed)
    return variants


def variant_overlaps(
    chrom: str, start: int, end: int, variants: Iterable[tuple[str, int, int, float | None]]
) -> bool:
    return any(_chrom(chrom) == vc and start < ve and vs < end for vc, vs, ve, _ in variants)


def _accession_keys(token: str) -> set[str]:
    """Exact token plus Ensembl/RefSeq-style identifier without trailing .version."""
    value = str(token).strip().lower()
    if not value:
        return set()
    keys = {value}
    base, dot, rest = value.rpartition(".")
    if dot and rest.isdigit() and base:
        keys.add(base)
    return keys


def _tokenize_subject(row: pd.Series) -> set[str]:
    text = " ".join(str(row.get(k, "")) for k in ("sseqid", "sallacc", "stitle")).lower()
    tokens: set[str] = set()
    for raw in re.findall(r"[a-z0-9_.-]+", text):
        tokens.update(_accession_keys(raw))
    return tokens


_BIOTYPE_LABELS = {
    "protein_coding": "protein-coding mRNA",
    "lncRNA": "lncRNA",
    "retained_intron": "retained intron",
    "nonsense_mediated_decay": "NMD transcript",
    "processed_transcript": "processed transcript",
    "protein_coding_CDS_not_defined": "protein-coding, CDS not defined",
    "miRNA": "miRNA",
    "snRNA": "snRNA",
    "snoRNA": "snoRNA",
    "rRNA": "rRNA",
    "rRNA_pseudogene": "rRNA pseudogene",
    "Mt_tRNA": "mitochondrial tRNA",
    "Mt_rRNA": "mitochondrial rRNA",
}


def _gencode_fields(title: str) -> dict[str, str]:
    """GENCODE header: transcript|gene|havana gene|havana transcript|name|gene name|length|biotype|."""
    parts = str(title).split("|")
    if len(parts) < 8:
        return {"gene": "", "transcript": "", "transcript_type": "", "transcript_length": ""}
    biotype = parts[7]
    return {
        "gene": parts[5],
        "transcript": parts[4] or parts[0],
        "transcript_type": _BIOTYPE_LABELS.get(biotype, biotype),
        "transcript_length": parts[6],
    }


def readable_blast_matches(hits: pd.DataFrame) -> pd.DataFrame:
    """One row per site and gene: the longest same-strand, non-self alignment."""
    columns = [
        "site", "matched_gene", "transcript", "transcript_type", "region",
        "identity_pct", "aligned_nt", "site_start", "site_end",
        "transcript_start", "transcript_end", "transcript_length",
        "n_transcripts", "evalue", "drops_site",
    ]
    if hits.empty or "stitle" not in hits.columns:
        return pd.DataFrame(columns=columns)
    rows = []
    for record in hits.to_dict("records"):
        if bool(record.get("is_self")) or not subject_same_orientation(record):
            continue
        fields = _gencode_fields(str(record.get("stitle", "")))
        sstart, send = int(record["sstart"]), int(record["send"])
        rows.append({
            "site": record["qseqid"],
            "matched_gene": fields["gene"],
            "transcript": fields["transcript"],
            "transcript_type": fields["transcript_type"],
            "region": str(record.get("match_region") or ""),
            "identity_pct": float(record["pident"]),
            "aligned_nt": int(record["length"]),
            "site_start": int(record["qstart"]),
            "site_end": int(record["qend"]),
            "transcript_start": min(sstart, send),
            "transcript_end": max(sstart, send),
            "transcript_length": fields["transcript_length"],
            "evalue": record["evalue"],
            "drops_site": bool(record.get("offtarget_risk")),
            "_key": fields["gene"] or fields["transcript"] or str(record.get("sseqid")),
        })
    if not rows:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(rows)
    counts = frame.groupby(["site", "_key"]).size().rename("n_transcripts")
    # Prefer a UTR hit that drops the site over a longer CDS-only alignment
    # for the same gene, so the readable row matches the filter decision.
    frame["_utr"] = frame["region"].str.contains("UTR", na=False)
    frame = frame.sort_values(
        ["drops_site", "_utr", "aligned_nt", "identity_pct"],
        ascending=[False, False, False, False],
    ).drop_duplicates(["site", "_key"])
    frame = frame.merge(counts.reset_index(), on=["site", "_key"])
    frame["drops_site"] = frame["drops_site"].map({True: "yes", False: "no"})
    frame = frame.sort_values(
        ["site", "drops_site", "aligned_nt"], ascending=[True, False, False]
    )
    return frame[columns]


def subject_same_orientation(row: pd.Series) -> bool:
    """True when the subject interval runs 5′→3′, same sequence as the sense site.

    blastn reports a minus-strand hit with sstart > send. That hit is the
    binding-element sequence itself on another transcript, not a place the
    element can pair.
    """
    return int(row["sstart"]) < int(row["send"])


def blast_hit_is_self(row: pd.Series, self_tokens: Iterable[str]) -> bool:
    subject_tokens = _tokenize_subject(row)
    return any(_accession_keys(token) & subject_tokens for token in self_tokens)


def blast_hit_is_risk(
    row: pd.Series,
    query_length: int,
    min_length: int = 20,
    min_identity: float = 0.0,
    min_coverage: float = 0.0,
) -> bool:
    """Flag a non-self hit as off-target risk.

    Default length gate is ≥ 20 nt (inclusive). Identity and coverage floors
    are optional and off by default so that abundant 3′UTR windows can be
    pruned early.
    """
    coverage = float(row["length"]) / max(1, query_length)
    return (
        int(row["length"]) >= min_length
        and float(row["pident"]) >= min_identity
        and coverage >= min_coverage
    )


def exclusion_reason(row: pd.Series) -> str:
    """Plain description of the UTR hit that removed a site."""
    fields = _gencode_fields(str(row.get("stitle", "")))
    gene = fields["gene"] or "unknown gene"
    transcript = fields["transcript"] or str(row.get("sseqid", ""))
    return (
        f"{row['match_region']} of {gene} ({transcript}), "
        f"{int(row['length'])} nt, {float(row['pident']):.1f}% identity, "
        f"site {int(row['qstart'])}-{int(row['qend'])}, "
        f"transcript {int(row['sstart'])}-{int(row['send'])}"
    )


def run_blast(
    candidates: pd.DataFrame,
    blast_db: str,
    logger: logging.Logger,
    self_tokens: Iterable[str],
    options: SelectOptions,
    blastn_exe: str | None = None,
    utr_index: dict | None = None,
) -> pd.DataFrame:
    exe = resolve_executable("blastn", blastn_exe)
    if not exe:
        raise RuntimeError("blastn_not_found")
    if utr_index is None:
        index_path = discover_utr_index()
        if index_path is None:
            raise RuntimeError(
                "GENCODE UTR annotation is required; run `bssf db init --gencode-v45-transcripts`."
            )
        utr_index = load_utr_index(index_path)
    query = "".join(
        f">{row['name']}\n{str(row['target_sequence']).upper().replace('U', 'T')}\n"
        for _, row in candidates.iterrows()
    )
    command = [
        exe, "-task", BLASTN_SHORT["task"], "-db", blast_db, "-query", "-",
        "-outfmt", "6 " + " ".join(BLAST_COLUMNS),
        "-max_target_seqs", str(BLASTN_SHORT["max_target_seqs"]),
    ]
    proc = subprocess.run(command, input=query, capture_output=True, text=True, timeout=300)
    log_command(logger, command, proc.returncode, proc.stdout, proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"blastn_exit_{proc.returncode}")
    rows = [line.split("\t") for line in proc.stdout.splitlines() if line.strip()]
    hits = pd.DataFrame(rows, columns=BLAST_COLUMNS)
    for column in ("pident", "length", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"):
        if column in hits:
            hits[column] = pd.to_numeric(hits[column])
    if not hits.empty:
        hits["is_self"] = hits.apply(lambda r: blast_hit_is_self(r, self_tokens), axis=1)
        hits["match_region"] = [
            utr_regions(utr_index or {}, str(title).split("|", 1)[0], int(start), int(end))
            if int(start) < int(end) else ""
            for title, start, end in zip(hits["stitle"], hits["sstart"], hits["send"])
        ]
        hits["offtarget_risk"] = hits.apply(
            lambda r: (
                not r["is_self"]
                and subject_same_orientation(r)
                and "UTR" in str(r["match_region"])
                and blast_hit_is_risk(
                    r, int(r["qlen"]), options.offtarget_min_length,
                    options.min_identity, options.min_coverage,
                )
            ),
            axis=1,
        )
    else:
        hits["is_self"] = pd.Series(dtype=bool)
        hits["offtarget_risk"] = pd.Series(dtype=bool)
    return hits


def _hash_target(target: Target) -> str:
    return hashlib.sha256(target.sequence.upper().encode()).hexdigest()


def _tool_version(exe: str | None, flag: str = "--version") -> str | None:
    if not exe:
        return None
    try:
        proc = subprocess.run([exe, flag], capture_output=True, text=True, timeout=10)
        return (proc.stdout or proc.stderr).strip().splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unknown"


def _package_versions() -> dict[str, str]:
    versions = {}
    for package in ("bssf", "pandas", "biopython", "typer", "requests"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def select(
    target: Target,
    options: SelectOptions | None = None,
    *,
    runs_dir: Path = Path("runs"),
    variant_vcf: str | Path | None = None,
    blast_db: str | None = None,
    sites: list[tuple[str, str]] | None = None,
    use_variants: bool = False,
    self_tokens: Iterable[str] = (),
    blastn_exe: str | None = None,
) -> Path:
    options = options or SelectOptions()
    now = datetime.now(timezone.utc)
    run_id = uuid4().hex[:8]
    run_dir = runs_dir / f"{now.strftime('%Y%m%dT%H%M%SZ')}_{run_id}"
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(parents=True)
    logger = run_logger(run_dir / "run.log")
    effective_config = load_config()
    manifest = {
        "run_id": run_id, "created_at": now.isoformat(), "status": "running",
        "model_version": MODEL_VERSION, "parameters": asdict(options),
        "input_sha256": _hash_target(target), "target": target.metadata(),
        "database": {
            "variant_vcf": str(variant_vcf) if variant_vcf else None,
            "blast_db": blast_db,
            "assembly": effective_config.get("assembly"),
            "variant_source": effective_config.get("variant_source"),
            "variant_release": effective_config.get("variant_release"),
            "transcriptome_release": effective_config.get("transcriptome_release"),
            "transcriptome_source": effective_config.get("transcriptome_source"),
            "transcriptome_assembly": effective_config.get("transcriptome_assembly"),
        },
        "command": sys.argv,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": _package_versions(),
        },
        "effective_config": effective_config,
        "tools": {},
        "stage_counts": {},
    }
    write_manifest(run_dir / "manifest.json", manifest)
    (inputs_dir / "target.fasta").write_text(f">{target.name}\n{target.sequence.upper()}\n")
    try:
        logger.info("stage=candidates")
        if sites is None:
            all_candidates = generate_windows(target, options.window_size, options.step)
        else:
            all_candidates = sites_on_utr(target, sites)
        manifest["stage_counts"]["generated"] = len(all_candidates)
        active = all_candidates["status"].eq("pending")

        logger.info("stage=blast")
        if not blast_db:
            raise RuntimeError("BLAST database is required; run `bssf db init`.")
        hits = pd.DataFrame(columns=BLAST_COLUMNS + ["is_self", "offtarget_risk"])
        if active.any():
            hits = run_blast(
                all_candidates.loc[active], blast_db, logger, self_tokens, options, blastn_exe
            )
            risky_hits = hits.loc[hits["offtarget_risk"]]
            for name, group in risky_hits.groupby("qseqid"):
                best = group.sort_values(["length", "pident"], ascending=False).iloc[0]
                all_candidates.loc[
                    all_candidates["name"].eq(name), ["status", "failure_reason"]
                ] = ["failed", exclusion_reason(best)]
            active = all_candidates["status"].eq("pending")
        manifest["stage_counts"]["after_blast"] = int(active.sum())
        hits.to_csv(run_dir / "blast_hits.tsv", sep="\t", index=False)
        readable_blast_matches(hits).to_csv(run_dir / "blast_matches.tsv", sep="\t", index=False)

        logger.info("stage=variants")
        if use_variants and variant_vcf:
            if target.chrom is None or target.start is None or target.end is None:
                logger.warning("variant filtering skipped: target has no genomic coordinates")
            else:
                variants = read_variant_vcf(
                    variant_vcf,
                    options.min_af,
                    chrom=str(target.chrom),
                    start=int(target.start),
                    end=int(target.end),
                )
                for idx, row in all_candidates.loc[active].iterrows():
                    if pd.isna(row["start"]) or pd.isna(row["end"]):
                        continue
                    if variant_overlaps(str(row["chrom"]), int(row["start"]), int(row["end"]), variants):
                        all_candidates.loc[idx, ["status", "failure_reason"]] = [
                            "failed", "variant_overlap",
                        ]
                active = all_candidates["status"].eq("pending")
        elif use_variants:
            logger.warning("variant filtering skipped: no VCF configured")
        else:
            logger.info("variant filtering off")
        manifest["stage_counts"]["after_variants"] = int(active.sum())

        all_candidates.loc[active, "status"] = "pass"
        passed = all_candidates["status"].eq("pass")

        for column in WINDOW_COLUMNS:
            if column not in all_candidates:
                all_candidates[column] = pd.NA
        passed_table = all_candidates.loc[passed, WINDOW_COLUMNS]
        all_candidates[WINDOW_COLUMNS].to_csv(run_dir / "all_binding_sites.tsv", sep="\t", index=False)
        passed_table.to_csv(run_dir / "candidates.tsv", sep="\t", index=False)
        manifest["status"] = "completed"
        manifest["candidate_count"] = int(passed.sum())
        manifest["stage_counts"]["passed"] = int(passed.sum())
        if not passed.any():
            logger.warning("pipeline completed with zero passing candidates")
        return run_dir
    except Exception as exc:
        logger.exception("pipeline failed")
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        raise
    finally:
        blastn = resolve_executable("blastn", blastn_exe)
        manifest["tools"] = {
            "blastn": _tool_version(blastn, "-version"),
        }
        write_manifest(run_dir / "manifest.json", manifest)
