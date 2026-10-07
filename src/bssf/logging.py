from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any


def run_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger(f"bssf.{path.parent.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.propagate = False
    return logger


def log_command(
    logger: logging.Logger,
    command: list[str],
    returncode: int,
    stdout: str,
    stderr: str,
) -> None:
    logger.info("command=%s", json.dumps(command))
    logger.info("returncode=%s", returncode)
    logger.info("stdout:\n%s", stdout)
    logger.info("stderr:\n%s", stderr)


def write_manifest(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, default=str) + "\n")
