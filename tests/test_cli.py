import json
from pathlib import Path

from typer.testing import CliRunner

from bsst.cli import _header_metadata, app

runner = CliRunner()


def test_legacy_fasta_header_metadata() -> None:
    assert _header_metadata(
        "target chrom=4 start=10 end=20 strand=-"
    ) == {"chrom": "4", "start": "10", "end": "20", "strand": "-"}


def test_help_check_requirements_and_config(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("BSST_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("PATH", str(tmp_path))
    assert runner.invoke(app, ["--help"]).exit_code == 0
    checked = runner.invoke(app, ["check_requirements"])
    assert checked.exit_code == 1
    assert "RNAup" in checked.stdout
    listed = runner.invoke(app, ["resources"])
    assert listed.exit_code == 0
    assert "GENCODE" in listed.stdout
    assert "NCBI dbSNP" in listed.stdout
    shown = runner.invoke(app, ["config", "show"])
    assert shown.exit_code == 0
    assert "runs_dir" in shown.stdout


def test_check_requirements_ok_with_required_binaries(tmp_path: Path, monkeypatch) -> None:
    for name in ("RNAup", "blastn", "makeblastdb"):
        exe = tmp_path / name
        exe.write_text("#!/bin/sh\nexit 0\n")
        exe.chmod(0o755)
    monkeypatch.setenv("BSST_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("PATH", str(tmp_path))
    checked = runner.invoke(app, ["check_requirements"])
    assert checked.exit_code == 0, checked.output
    assert "missing/optional" in checked.stdout


def test_filter_help_requires_gene() -> None:
    result = runner.invoke(app, ["filter", "--help"])
    assert result.exit_code == 0
    assert "--gene" in result.stdout
    assert "--sites" in result.stdout
    assert "--rnaup" in result.stdout
    missing = runner.invoke(app, ["filter", "--utr", "examples/LETM1.fasta"])
    assert missing.exit_code != 0
