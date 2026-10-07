# Local scientific data (not Git-tracked)

Pinned downloads named in `src/bsst/resources.py`. Fetch and verify with
`bsst db init --dbsnp-common-all --gencode-v45-transcripts` and
`bsst check_requirements`. Details: [docs/resources.md](../docs/resources.md).

| File | Producer |
|------|----------|
| `dbSNP_b151_GRCh38p7_common_all_20180418.vcf.gz` (~1.5 GB) | NCBI dbSNP b151, GRCh38.p7, `common_all` 20180418 |
| `gencode.v45.transcripts.fa` (~454 MB) + `gencode_v45_transcripts_db.*` | GENCODE 45 / Ensembl 111 / GRCh38.p14 / CHR |
| `gencode.v45.annotation.gff3.gz` → `gencode.v45.utr_on_transcript.tsv` | same GENCODE release; 5′/3′UTR on spliced transcripts |

dbSNP b151 is GRCh38.**p7**. GENCODE 45 / Ensembl 111 are GRCh38.**p14**.
Each run `manifest.json` records those fields separately.
