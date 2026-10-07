# Binding Site Safety Filter (`bssf`) guide

`bssf` is a **safety filter** for Hepha antisense sites on a transcript
3′UTR: drop same-strand UTR offtargets (BLAST), optionally drop common-SNP
overlaps (`--variants`). A pass is not upregulation.

Config lives under `~/.bssf` (override with `BSSF_HOME`).

```text
3′UTR (--utr or Ensembl 111) → sites (--sites or 40 nt windows)
  → BLAST GENCODE 45 (same-strand, non-self, ≥20 nt on 5′/3′UTR → drop)
  → optional dbSNP COMMON via tabix (--variants)
  → candidates.tsv / all_binding_sites.tsv
```

`--gene` is required so self hits are ignored. Prefer `--sites`; omitting
it blasts every window in one 300 s `blastn` and will time out on long UTRs.
`--utr` is the reproducible path; omitting it uses Ensembl archive 111’s
canonical coding transcript only.

## BLAST

`blastn -task blastn-short` vs GENCODE 45 CHR transcripts (GRCh38.p14,
252930 seqs). Queries are DNA (`T`). A site fails when a same-strand,
non-self hit ≥ `--offtarget-min-length` (default 20) overlaps a 5′UTR or
3′UTR on the spliced transcript. CDS-only hits, reverse-strand hits, and
transcripts with no UTR annotation stay. Reverse-complement subject spans
(`sstart > send`) are ignored. Do not pass `-parse_seqids` to
`makeblastdb` (GENCODE headers contain `|`).

Read `blast_matches.tsv`: one row per site/gene; `drops_site=yes` only for
UTR hits that also meet the length cutoff.

## Variants

Off unless `--variants`. NCBI dbSNP **b157** on **GRCh38.p14**
(`GCF_000001405.40`). `bssf db init --dbsnp-common-all` downloads only the
**~3 MB tabix index**; `--variants` runs `tabix` against the HTTPS VCF for
the target UTR interval, keeps `INFO/COMMON` sites, and remaps RefSeq
contigs (`NC_000001.11` → `1`, …). Needs UTR `chrom` + 0-based BED
`start`/`end`. Optional `--min-af`. Requires `tabix` on `PATH` and network
access to NCBI (a local HTTP proxy helps).

(`p` in GRCh38.**p14** is the GRC **patch** level of the same major assembly.)

## Outputs (`runs/<UTC>_<id>/`)

| File | Content |
|------|---------|
| `candidates.tsv` | `status=pass` |
| `all_binding_sites.tsv` | every site + `failure_reason` |
| `blast_matches.tsv` | offtarget summary |
| `blast_hits.tsv` | raw BLAST |
| `run.log` / `manifest.json` | commands, versions, DB pins |

Typical `failure_reason`: `not_in_utr`, `ambiguous`, `variant_overlap`, or
a UTR match like `3'UTR of GENE (tx), 21 nt, 95% identity, …`.

## Pinned resources

Identities live in `src/bssf/resources.py` (`bssf resources`,
`bssf check_requirements`). BLAST, UTR index, Ensembl fetch, and variants
all pin **GRCh38.p14**.

| Role | Pin |
|------|-----|
| Variants | dbSNP b157 `GCF_000001405.40` via tabix (COMMON), GRCh38.p14 |
| BLAST | GENCODE 45 CHR transcripts, GRCh38.p14 |
| UTR index | same GFF3 → `gencode.v45.utr_on_transcript.tsv` |
| Fetch | `https://e111.rest.ensembl.org` |
| Examples | `examples/LETM1.*`, `examples/NSD2.*` |

```bash
bssf db init --dbsnp-common-all --gencode-v45-transcripts
```

## References

- Camacho et al. (2009). BLAST+: architecture and applications. *BMC
  Bioinformatics* **10**, 421. https://doi.org/10.1186/1471-2105-10-421
- Frankish et al. (2023). GENCODE: reference annotation for the human and
  mouse genomes in 2023. *Nucleic Acids Research* **51**, D942–D949.
  https://doi.org/10.1093/nar/gkac1071
- Harrison et al. (2024). Ensembl 2024. *Nucleic Acids Research* **52**,
  D891–D899. https://doi.org/10.1093/nar/gkad1049
- Sherry et al. (2001). dbSNP: the NCBI database of genetic variation.
  *Nucleic Acids Research* **29**, 308–311.
  https://doi.org/10.1093/nar/29.1.308
- Li (2011). Tabix: fast retrieval of sequence features from generic
  TAB-delimited files. *Bioinformatics* **27**, 718–719.
  https://doi.org/10.1093/bioinformatics/btq671
