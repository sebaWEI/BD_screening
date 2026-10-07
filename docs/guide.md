# Binding Site Safety Filter (`bssf`) Guide

`bssf` is a **sequence-level safety and robustness filter** for HEPHA antisense binding sites in transcript 3′UTRs. It removes candidate sites with potentially problematic same-strand UTR off-targets identified by BLAST and, optionally, sites overlapping common dbSNP variants.

**A site passing `bssf` is sequence-level eligible, not functionally optimal.** Functional prioritization is performed by our model-guided BD design framework.

Configuration is stored under `~/.bssf` (override with `BSSF_HOME`).

```text
3′UTR (--utr or Ensembl 111)
    → candidate sites (--sites or 40 nt windows)
    → BLAST against GENCODE 45
       (same-strand, non-self, ≥20 nt UTR hit → drop)
    → optional dbSNP COMMON filtering via tabix (--variants)
    → candidates.tsv / all_binding_sites.tsv
```

`--gene` is required so that self hits are ignored. Prefer `--sites` when possible; omitting it causes `bssf` to BLAST every 40 nt window in a single 300 s `blastn` run, which may time out for long UTRs.

`--utr` is the recommended reproducible input path. If omitted, `bssf` uses the canonical coding transcript from the archived Ensembl 111 endpoint as a convenience fallback. This fallback does not change the pinned BLAST or variant resources.

## BLAST

`bssf` uses `blastn -task blastn-short` against the **GENCODE 45 CHR transcript set** on **GRCh38.p14** (252,930 sequences). Candidate sites are represented as DNA sequences (`T` rather than `U`).

A site fails when a same-strand, non-self BLAST hit of at least `--offtarget-min-length` (default: 20 nt) overlaps a 5′UTR or 3′UTR in the spliced transcript.

CDS-only hits, reverse-strand hits, and transcripts without UTR annotation do not cause a candidate to be dropped. Reverse-complement subject spans (`sstart > send`) are ignored.

Do not pass `-parse_seqids` to `makeblastdb`, because GENCODE FASTA headers contain `|` characters.

`blast_matches.tsv` contains one row per site/gene pair. `drops_site=yes` is assigned only when a hit overlaps an annotated UTR and satisfies the minimum length threshold.

## Variants

Variant filtering is disabled by default and can be enabled with `--variants`.

`bssf` uses **NCBI dbSNP b157** on **GRCh38.p14** (`GCF_000001405.40`). The initialization command downloads only the approximately 3 MB Tabix index. When `--variants` is enabled, `bssf` queries the remote VCF through `tabix` for the target UTR interval, retains records with `INFO/COMMON`, and remaps RefSeq chromosome identifiers (`NC_000001.11` → `1`, etc.).

The UTR interval must provide `chrom` together with 0-based BED `start` and `end` coordinates. An optional `--min-af` threshold can be applied.

`tabix` must be available on `PATH`, and network access to the NCBI VCF endpoint is required. Set `https_proxy` if needed; otherwise `bssf` auto-detects a local HTTP proxy on `127.0.0.1:7897` / `7890` / `1087` / `8080`.

The `p` in `GRCh38.p14` denotes the GRC **patch level** of the GRCh38 reference assembly.

## Outputs

Results are stored under `runs/<UTC>_<id>/`.

| File | Content |
|------|---------|
| `candidates.tsv` | Candidates with `status=pass` |
| `all_binding_sites.tsv` | All evaluated sites with `failure_reason` |
| `blast_matches.tsv` | BLAST off-target summary |
| `blast_hits.tsv` | Raw BLAST hits |
| `run.log` / `manifest.json` | Commands, software versions, parameters, and database pins |

Typical `failure_reason` values include `not_in_utr`, `ambiguous`, `variant_overlap`, or a UTR match such as `3'UTR of GENE (tx), 21 nt, 95% identity, …`.

## Pinned Resources

Resource identities are defined in `src/bssf/resources.py` and can be inspected with:

```bash
bssf resources
bssf check_requirements
```

BLAST, UTR indexing, and variant resources are pinned to **GRCh38.p14**. The Ensembl fallback uses the archived Ensembl 111 endpoint.

| Role | Pin |
|------|-----|
| Variants | dbSNP b157, `GCF_000001405.40`, COMMON subset, GRCh38.p14 |
| BLAST | GENCODE 45 CHR transcripts, GRCh38.p14 |
| UTR index | GENCODE 45 GFF3 → `gencode.v45.utr_on_transcript.tsv` |
| Fetch fallback | Ensembl archive 111 |
| Examples | `examples/LETM1.*`, `examples/NSD2.*` |

Initialize the local resources with:

```bash
bssf db init --dbsnp-common-all --gencode-v45-transcripts
```

## References

- Camacho et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* **10**, 421. https://doi.org/10.1186/1471-2105-10-421
- Frankish et al. (2023). GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Research* **51**, D942–D949. https://doi.org/10.1093/nar/gkac1071
- Harrison et al. (2024). Ensembl 2024. *Nucleic Acids Research* **52**, D891–D899. https://doi.org/10.1093/nar/gkad1049
- Sherry et al. (2001). dbSNP: the NCBI database of genetic variation. *Nucleic Acids Research* **29**, 308–311. https://doi.org/10.1093/nar/29.1.308
- Li (2011). Tabix: fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics* **27**, 718–719. https://doi.org/10.1093/bioinformatics/btq671