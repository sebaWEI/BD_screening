from pathlib import Path

import pandas as pd

from bssf.models import Target
from bssf.pipeline import generate_windows, readable_blast_matches, select, sites_on_utr


def _empty_blast(*_args, **_kwargs) -> pd.DataFrame:
    return pd.DataFrame({"qseqid": [], "offtarget_risk": []})


def test_unique_site_is_placed_on_utr() -> None:
    target = Target("AACCGGTTAA", "utr", "4", 100, 110, "+")
    placed = sites_on_utr(target, [("site", "CCGG")])
    assert placed.loc[0, "status"] == "pending"
    assert placed.loc[0, "utr_start"] == 2
    assert placed.loc[0, "start"] == 102
    assert placed.loc[0, "end"] == 106


def test_missing_and_repeated_sites_are_not_pending() -> None:
    target = Target("AAAA", "utr")
    placed = sites_on_utr(target, [("gone", "CC"), ("repeat", "AA")])
    assert placed.loc[0, "failure_reason"] == "not_in_utr"
    assert placed.loc[1, "failure_reason"] == "ambiguous"
    assert placed["status"].eq("pending").sum() == 0


def test_omitted_sites_slide_and_variant_stage_stays_off(tmp_path: Path, monkeypatch) -> None:
    calls = {"variants": 0}

    monkeypatch.setattr("bssf.pipeline.run_blast", _empty_blast)
    monkeypatch.setattr(
        "bssf.pipeline.read_variant_vcf",
        lambda *_args, **_kwargs: calls.__setitem__("variants", 1),
    )
    target = Target("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA", "utr", gene="GENE")
    run_dir = select(
        target,
        runs_dir=tmp_path,
        blast_db="unused",
        self_tokens=["GENE"],
    )
    table = pd.read_csv(run_dir / "candidates.tsv", sep="\t")
    assert len(table) == len(generate_windows(target, 40))
    assert set(table["status"]) == {"pass"}
    assert calls == {"variants": 0}
    log = (run_dir / "run.log").read_text()
    assert "variant filtering off" in log
    assert "stage=rnaup" not in log


def test_blast_drop_keeps_only_passing_sites(tmp_path: Path, monkeypatch) -> None:
    def blast(candidates, *_args, **_kwargs):
        return pd.DataFrame([{
            "qseqid": "bad", "offtarget_risk": True, "length": 20, "pident": 100.0,
            "qstart": 1, "qend": 20, "sstart": 10, "send": 29,
            "stitle": "ENST1|ENSG1|h1|h2|OTHER-201|OTHER|100|protein_coding|",
            "sseqid": "ENST1", "match_region": "3'UTR", "evalue": 0.5,
        }])

    monkeypatch.setattr("bssf.pipeline.run_blast", blast)
    target = Target("AACCGGTT", "utr", gene="GENE")
    run_dir = select(
        target,
        runs_dir=tmp_path,
        blast_db="unused",
        sites=[("keep", "AACC"), ("bad", "GGTT"), ("gone", "TTTT")],
        self_tokens=["GENE"],
    )
    all_rows = pd.read_csv(run_dir / "all_binding_sites.tsv", sep="\t")
    reasons = dict(zip(all_rows["name"], all_rows["failure_reason"].fillna(""), strict=True))
    assert reasons["gone"] == "not_in_utr"
    assert reasons["bad"].startswith("3'UTR of OTHER")
    passed = pd.read_csv(run_dir / "candidates.tsv", sep="\t")
    assert passed["name"].tolist() == ["keep"]
    assert set(passed["status"]) == {"pass"}


def test_readable_blast_matches_collapse_to_one_gene() -> None:
    title = (
        "ENST000001.1|ENSG000001.1|OTTHUMG1|OTTHUMT1|"
        "OXCT1-201|OXCT1|3000|protein_coding|"
    )
    other = title.replace("OXCT1-201", "OXCT1-202").replace("ENST000001.1", "ENST000002.1")
    hits = pd.DataFrame([
        {
            "qseqid": "site", "stitle": title, "sseqid": "a", "pident": 92.3,
            "length": 26, "qstart": 7, "qend": 32, "sstart": 80, "send": 100,
            "evalue": 0.13, "is_self": False, "offtarget_risk": True,
        },
        {
            "qseqid": "site", "stitle": other, "sseqid": "b", "pident": 100.0,
            "length": 20, "qstart": 1, "qend": 20, "sstart": 1, "send": 20,
            "evalue": 1.0, "is_self": False, "offtarget_risk": False,
        },
        {
            "qseqid": "site", "stitle": title, "sseqid": "self", "pident": 100.0,
            "length": 40, "qstart": 1, "qend": 40, "sstart": 1, "send": 40,
            "evalue": 1e-10, "is_self": True, "offtarget_risk": False,
        },
        {
            "qseqid": "site", "stitle": title.replace("OXCT1", "RBM44"), "sseqid": "rev",
            "pident": 100.0, "length": 40, "qstart": 1, "qend": 40,
            "sstart": 200, "send": 161, "evalue": 1e-8,
            "is_self": False, "offtarget_risk": True,
        },
    ])
    table = readable_blast_matches(hits)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["matched_gene"] == "OXCT1"
    assert row["transcript_type"] == "protein-coding mRNA"
    assert row["aligned_nt"] == 26
    assert row["transcript_start"] == 80
    assert row["transcript_end"] == 100
    assert row["n_transcripts"] == 2
    assert row["drops_site"] == "yes"


def test_readable_blast_matches_prefers_utr_drop_over_longer_cds() -> None:
    base = (
        "ENST000001.1|ENSG000001.1|OTTHUMG1|OTTHUMT1|"
        "OTHER-201|OTHER|3000|protein_coding|"
    )
    hits = pd.DataFrame([
        {
            "qseqid": "site", "stitle": base, "sseqid": "cds", "pident": 100.0,
            "length": 30, "qstart": 1, "qend": 30, "sstart": 500, "send": 529,
            "evalue": 0.01, "is_self": False, "offtarget_risk": False,
            "match_region": "",
        },
        {
            "qseqid": "site", "stitle": base.replace("OTHER-201", "OTHER-202"),
            "sseqid": "utr", "pident": 95.0, "length": 22, "qstart": 5, "qend": 26,
            "sstart": 10, "send": 31, "evalue": 0.2, "is_self": False,
            "offtarget_risk": True, "match_region": "3'UTR",
        },
    ])
    table = readable_blast_matches(hits)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["region"] == "3'UTR"
    assert row["aligned_nt"] == 22
    assert row["drops_site"] == "yes"
