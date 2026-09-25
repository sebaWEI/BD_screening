from __future__ import annotations

from time import sleep

import requests

from .models import Target
from .resources import ENSEMBL_REST


class EnsemblArchiveError(RuntimeError):
    """The pinned Ensembl REST archive is missing or is not release 111."""


def _ensembl_get(path: str, *, content_type: str, params: dict[str, str] | None = None) -> requests.Response:
    """Ensembl archive REST requires content-type as a query parameter."""
    query = {"content-type": content_type}
    if params:
        query.update(params)
    last_exc: requests.RequestException | None = None
    attempts = 4
    for attempt in range(attempts):
        try:
            response = requests.get(
                f"{ENSEMBL_REST['base_url']}{path}",
                params=query,
                timeout=60,
            )
            response.raise_for_status()
            return response
        except requests.HTTPError as exc:
            last_exc = exc
            status = exc.response.status_code if exc.response is not None else 0
            if status and status < 500 and status != 429:
                break
        except requests.RequestException as exc:
            last_exc = exc
        if attempt + 1 < attempts:
            sleep(0.5 * (2 ** attempt))
    raise EnsemblArchiveError(
        f"Ensembl REST archive {ENSEMBL_REST['base_url']} request failed for {path}. "
        "Pass `--utr` with a stored 3'UTR. "
        "bsst does not fall back to rest.ensembl.org."
    ) from last_exc


_cached_release: str | None = None


def assert_ensembl_111() -> str:
    """Refuse to run against a live or wrong-release REST server."""
    global _cached_release
    if _cached_release is not None:
        return _cached_release
    payload = _ensembl_get("/info/software", content_type="application/json").json()
    release = str(payload.get("release", ""))
    expected = str(ENSEMBL_REST["ensembl_version"])
    if release != expected:
        raise EnsemblArchiveError(
            f"{ENSEMBL_REST['base_url']} reported Ensembl release {release}, "
            f"expected {expected} (GENCODE {ENSEMBL_REST['gencode_release']}). "
            "Pass `--utr` with a stored 3'UTR."
        )
    _cached_release = release
    return release


def fetch_region(
    chrom: str,
    start_1: int,
    end_1: int,
    strand: int,
    species: str | None = None,
) -> str:
    """Fetch a 1-based inclusive genomic interval from Ensembl 111.

    ``strand`` is ``1`` (plus) or ``-1`` (minus, reverse-complemented by Ensembl).
    """
    assert_ensembl_111()
    species = species or ENSEMBL_REST["default_species"]
    text = _ensembl_get(
        f"/sequence/region/{species}/{chrom}:{start_1}..{end_1}:{strand}",
        content_type="text/plain",
    ).text.strip().upper().replace("U", "T")
    if not text:
        raise EnsemblArchiveError(f"empty sequence for {chrom}:{start_1}..{end_1}:{strand}")
    return text


def fetch_bed_sequence(chrom: str, start: int, end: int, strand: str) -> str:
    """Fetch a 0-based half-open BED interval as transcript-oriented DNA."""
    if strand not in {"+", "-"}:
        raise ValueError("strand must be + or -")
    strand_i = 1 if strand == "+" else -1
    return fetch_region(chrom, start + 1, end, strand_i)


def fetch_gene_utr(
    gene: str, species: str | None = None,
) -> Target:
    """Fetch the canonical coding transcript 3'UTR from Ensembl 111 REST."""
    release = assert_ensembl_111()
    species = species or ENSEMBL_REST["default_species"]
    data = _ensembl_get(
        f"/lookup/symbol/{species}/{gene}",
        content_type="application/json",
        params={"expand": "1"},
    ).json()
    transcripts = data.get("Transcript", [])
    canonical = next((t for t in transcripts if t.get("is_canonical") == 1), None)
    if canonical is None and transcripts:
        canonical = max(transcripts, key=lambda t: t["end"] - t["start"])
    if not canonical or "Translation" not in canonical:
        raise ValueError(f"{gene} has no canonical coding transcript")

    strand_i = int(canonical["strand"])
    if strand_i == 1:
        start_1, end_1 = canonical["Translation"]["end"] + 1, canonical["end"]
    else:
        start_1, end_1 = canonical["start"], canonical["Translation"]["start"] - 1
    if start_1 > end_1:
        raise ValueError(f"{gene} has no non-empty 3'UTR")

    chrom = str(data["seq_region_name"])
    sequence = fetch_region(chrom, start_1, end_1, strand_i, species)
    return Target(
        sequence=sequence,
        name=f"{gene}_{canonical['id']}",
        chrom=chrom,
        start=start_1 - 1,
        end=end_1,
        strand="+" if strand_i == 1 else "-",
        gene=gene,
        transcript_id=canonical["id"],
        assembly=data.get("assembly_name"),
        annotation_release=release,
        source=(
            f"Ensembl REST archive {ENSEMBL_REST['base_url']} "
            f"(Ensembl {release} / GENCODE {ENSEMBL_REST['gencode_release']}, {species})"
        ),
    )
