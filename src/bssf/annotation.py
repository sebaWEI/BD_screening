"""Map a spliced-transcript BLAST interval onto GENCODE 5′UTR and 3′UTR."""

from __future__ import annotations

import gzip
from pathlib import Path


def _attrs(field: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in field.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def build_utr_index(gff3: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Project GENCODE UTR features onto spliced transcript coordinates."""
    exons: dict[str, list[tuple[int, int, str, int]]] = {}
    utrs: dict[str, list[tuple[int, int, str]]] = {}
    with _open_text(gff3) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature, strand = parts[2], parts[6]
            if feature not in {"exon", "five_prime_UTR", "three_prime_UTR"}:
                continue
            attrs = _attrs(parts[8])
            transcript_id = attrs.get("transcript_id") or ""
            if not transcript_id:
                continue
            start, end = int(parts[3]), int(parts[4])
            if feature == "exon":
                exons.setdefault(transcript_id, []).append(
                    (start, end, strand, int(attrs.get("exon_number") or 0))
                )
            else:
                region = "5'UTR" if feature == "five_prime_UTR" else "3'UTR"
                utrs.setdefault(transcript_id, []).append((start, end, region))
    index: dict[str, list[tuple[int, int, str]]] = {}
    for transcript_id, features in utrs.items():
        mapped = _project(exons.get(transcript_id, []), features)
        if mapped:
            index[transcript_id] = _merge(mapped)
    return index


def _project(
    exons: list[tuple[int, int, str, int]],
    utrs: list[tuple[int, int, str]],
) -> list[tuple[int, int, str]]:
    if not exons:
        return []
    ordered = sorted(exons, key=lambda item: item[3] or item[0])
    layout: list[tuple[int, int, str, int]] = []
    cursor = 1
    for genomic_start, genomic_end, strand, _number in ordered:
        length = genomic_end - genomic_start + 1
        layout.append((genomic_start, genomic_end, strand, cursor))
        cursor += length
    projected: list[tuple[int, int, str]] = []
    for utr_start, utr_end, region in utrs:
        for genomic_start, genomic_end, strand, spliced_start in layout:
            overlap_start = max(utr_start, genomic_start)
            overlap_end = min(utr_end, genomic_end)
            if overlap_start > overlap_end:
                continue
            if strand == "-":
                spliced_a = spliced_start + (genomic_end - overlap_end)
                spliced_b = spliced_start + (genomic_end - overlap_start)
            else:
                spliced_a = spliced_start + (overlap_start - genomic_start)
                spliced_b = spliced_start + (overlap_end - genomic_start)
            projected.append((min(spliced_a, spliced_b), max(spliced_a, spliced_b), region))
    return projected


def _merge(intervals: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    merged: list[tuple[int, int, str]] = []
    for start, end, region in sorted(intervals):
        if merged and merged[-1][2] == region and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end), region)
        else:
            merged.append((start, end, region))
    return merged


def write_utr_index(gff3: Path, destination: Path) -> Path:
    index = build_utr_index(gff3)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = ["transcript_id\tstart\tend\tregion\n"]
    for transcript_id in sorted(index):
        for start, end, region in index[transcript_id]:
            lines.append(f"{transcript_id}\t{start}\t{end}\t{region}\n")
    destination.write_text("".join(lines), encoding="utf-8")
    return destination


def load_utr_index(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    index: dict[str, list[tuple[int, int, str]]] = {}
    with path.open(encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            transcript_id, start, end, region = line.rstrip("\n").split("\t")
            index.setdefault(transcript_id, []).append((int(start), int(end), region))
    return index


def utr_regions(
    index: dict[str, list[tuple[int, int, str]]],
    transcript_id: str,
    start: int,
    end: int,
) -> str:
    """Return ``5'UTR``, ``3'UTR``, ``5'UTR+3'UTR``, or ``""`` when the interval misses both."""
    intervals = index.get(transcript_id)
    if intervals is None and "." in transcript_id:
        base = transcript_id.rsplit(".", 1)[0]
        intervals = next(
            (value for key, value in index.items() if key == base or key.startswith(base + ".")),
            None,
        )
    if not intervals:
        return ""
    found: list[str] = []
    for interval_start, interval_end, region in intervals:
        if start <= interval_end and end >= interval_start and region not in found:
            found.append(region)
    return "+".join(found)
