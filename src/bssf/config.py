from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


def config_dir() -> Path:
    return Path(os.environ.get("BSSF_HOME", Path.home() / ".bssf"))


def config_path() -> Path:
    return config_dir() / "config.json"


def default_config() -> dict[str, Any]:
    return {
        "db_dir": str(config_dir() / "db"),
        "variant_vcf": None,
        "blast_db": None,
        "blastn_exe": None,
        "assembly": None,
        "variant_source": None,
        "variant_release": None,
        "transcriptome_release": None,
        "transcriptome_source": None,
        "transcriptome_assembly": None,
        "runs_dir": "runs",
    }


def load_config() -> dict[str, Any]:
    cfg = default_config()
    path = config_path()
    if path.exists():
        cfg.update(json.loads(path.read_text()))
    return cfg


def save_config(config: dict[str, Any]) -> Path:
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2) + "\n")
    return path
