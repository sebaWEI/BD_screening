# Scientific limitations

This tool filters sense sites on the **3′UTR**. That is how Hepha is
specified. Passing BLAST (and optionally the variant gate) is a
**safety screen**, not proof of up-regulation. A default run does not
score, rank, or predict translational activity.

## What the filters omit

BLAST always uses the pinned transcriptome and the GENCODE 45 UTR index.
Without `data/gencode.v45.utr_on_transcript.tsv`, `filter` refuses to run.
A site is dropped only when a same-strand hit of at least
`--offtarget-min-length` nt (default 20) overlaps a 5′UTR or 3′UTR on the
spliced transcript. Hits on the reverse strand, CDS-only hits, and
transcripts with no UTR annotation are kept. The variant filter runs only
with `--variants` (and optional `--min-af`).

Neither filter predicts whether a site will upregulate translation. They
only ask whether the antisense sequence has a same-strand UTR offtarget
(and, optionally, a common SNP inside the window).

**`--utr`** is the reproducible path. Omitting it always takes the
Ensembl 111 archive’s canonical coding transcript, which may not be the
isoform you want, and the archive will eventually retire.

## Pinned data (not interchangeable)

| Resource | What it actually is |
|----------|---------------------|
| BLAST subject | GENCODE 45 / Ensembl 111 / GRCh38.**p14** / all CHR transcripts |
| UTR coordinates | same release, `gencode.v45.annotation.gff3.gz`, projected to `gencode.v45.utr_on_transcript.tsv` |
| Variant VCF | NCBI dbSNP b151 `common_all` / GRCh38.**p7** |
| omitted `--utr` | Ensembl REST **archive 111** only (`e111.rest.ensembl.org`) |

Catalog and verification: [resources.md](resources.md).

## After filtering

Review sites that passed, then test in the lab: negative controls, dose
response, RNA and protein readouts, and replication. Do not pool pre-filter
result tables with current pass/fail runs.

Tool and database citations: [methods.md](methods.md#references).
