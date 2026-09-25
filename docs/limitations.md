# Scientific limitations

This tool filters sense sites on the **3′UTR**. That is how Hepha is
specified. Passing BLAST is not proof of up-regulation. A default run
does not score or rank sites. `--rnaup` adds a thermodynamic hypothesis
for sites that already passed.

## What `--rnaup` still omits

RNAup is off unless you pass `--rnaup`. It is a thermodynamic RNA–RNA
model, and it does not know about:

- protein occupancy or RBP maps
- cellular compartment, modifications, or degradation
- isoform abundance or construct context
- genome (DNA) off-targets — BLAST here is against **transcripts**

BLAST always uses the pinned transcriptome and the GENCODE 45 UTR index.
Without `data/gencode.v45.utr_on_transcript.tsv`, `filter` refuses to run.
A site is dropped only when a same-strand hit of at least
`--offtarget-min-length` nt (default 20) overlaps a 5′UTR or 3′UTR on the
spliced transcript. Hits on the reverse strand, CDS-only hits, and
transcripts with no UTR annotation are kept. The variant filter runs only
with `--variants` (and optional `--min-af`).

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
response, RNA and protein readouts, and replication. A ΔG rank from
`--rnaup` is not a substitute for that. Old result files ranked every
window by RNAup; do not pool them with current pass/fail tables.

Tool and database citations: [methods.md](methods.md#references).
