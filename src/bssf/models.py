from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


MODEL_VERSION = "utr3-blast-variant-v1"


@dataclass(frozen=True)
class Target:
    sequence: str
    name: str
    chrom: str | None = None
    start: int | None = None
    end: int | None = None
    strand: str = "+"
    gene: str | None = None
    transcript_id: str | None = None
    assembly: str | None = None
    annotation_release: str | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        sequence = self.sequence.upper().replace("U", "T")
        if not sequence or set(sequence) - set("ACGTN"):
            raise ValueError("target sequence must contain only A/C/G/T/U/N")
        if self.strand not in {"+", "-"}:
            raise ValueError("strand must be + or -")
        if (self.start is None) != (self.end is None):
            raise ValueError("start and end must be supplied together")
        if self.start is not None and (
            self.start < 0 or self.end is None or self.end <= self.start
        ):
            raise ValueError("coordinates must be a valid 0-based half-open interval")
        object.__setattr__(self, "sequence", sequence)

    def metadata(self) -> dict[str, Any]:
        data = asdict(self)
        data.pop("sequence")
        data["length"] = len(self.sequence)
        data["mode"] = "3UTR"
        return data


@dataclass(frozen=True)
class SelectOptions:
    window_size: int = 40
    step: int = 1
    min_gc: float = 0.3
    max_gc: float = 0.7
    max_homopolymer: int = 4
    min_af: float | None = None
    offtarget_min_length: int = 20
    min_identity: float = 0.0
    min_coverage: float = 0.0

    def __post_init__(self) -> None:
        if self.window_size <= 0 or self.step <= 0:
            raise ValueError("window_size/step must be positive")
        if not 0 <= self.min_gc <= self.max_gc <= 1:
            raise ValueError("GC thresholds must satisfy 0 <= min_gc <= max_gc <= 1")
        if self.min_af is not None and not 0 <= self.min_af <= 1:
            raise ValueError("min_af must be between 0 and 1")
        if self.offtarget_min_length <= 0:
            raise ValueError("offtarget_min_length must be positive")
        if not 0 <= self.min_identity <= 100 or not 0 <= self.min_coverage <= 1:
            raise ValueError("BLAST identity/coverage thresholds are out of range")
