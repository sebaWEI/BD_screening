# Scientific resources

Pinned names live in `src/bsst/resources.py`. Print the same catalog with
`bsst resources`. `bsst check_requirements` checks local files and binaries
against those identities.

The variant VCF and the GENCODE FASTA are **different GRCh38 patches**.
The pipeline records each producer string separately.


| Role           | File / endpoint                         | Version                                | Assembly    |
| -------------- | --------------------------------------- | -------------------------------------- | ----------- |
| Variant VCF    | NCBI dbSNP `common_all_20180418.vcf.gz` | b151, 20180418                         | GRCh38.p7   |
| BLAST subject  | GENCODE `gencode.v45.transcripts.fa.gz` | Release 45 = Ensembl 111, 2024-01, CHR | GRCh38.p14  |
| 3′UTR fetch    | `https://e111.rest.ensembl.org`         | Ensembl 111                            | GRCh38.p14  |
| Off-target     | NCBI BLAST+ `blastn -task blastn-short` | local `blastn -version`                | —           |
| Binding energy | ViennaRNA `RNAup`                       | local `RNAup --version`                | —           |
| Example LETM1  | `examples/LETM1.fasta`                  | ENST00000302787 (v45 LETM1-201 `.3`)   | chrom=4 BED |
| Example NSD2   | `examples/NSD2.FASTA`                   | ENST00000508803 (v45 NSD2-218 `.6`)    | chrom=4 BED |


`--utr` is the durable path. Omitting it talks only to the Ensembl 111
archive. The bundled FASTAs match that archive’s canonical 3′UTRs for
those BED intervals; minus-strand LETM1 windows equal the reverse
complement of the plus-strand genome slice.

## Checksums

`bsst db init --dbsnp-common-all --gencode-v45-transcripts` checks SHA-256
of the downloaded bytes (and the uncompressed GENCODE FASTA) against
`src/bsst/resources.py`. Those digests match the producer MD5 files (NCBI
`*.vcf.gz.md5`, GENCODE `MD5SUMS`). A mismatch aborts; `--force`
re-downloads. Custom `--variant-url` / `--transcriptome-url` need an
explicit `--*-sha256` if you want a checksum.

```bash
uv run bsst db init --dbsnp-common-all --gencode-v45-transcripts
```

Already on disk:

```bash
uv run bsst db init \
  --variant-vcf data/dbSNP_b151_GRCh38p7_common_all_20180418.vcf.gz \
  --assembly GRCh38.p7 \
  --variant-source "NCBI dbSNP" \
  --variant-release "b151 fileDate=20180418 common_all" \
  --blast-db data/gencode_v45_transcripts_db \
  --transcriptome-source GENCODE \
  --transcriptome-release "GENCODE 45 / Ensembl 111 / 2024-01 / CHR transcripts" \
  --transcriptome-assembly GRCh38.p14
```



## Variant VCF (NCBI dbSNP b151, GRCh38.p7)

[https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz)

Local: `data/dbSNP_b151_GRCh38p7_common_all_20180418.vcf.gz` (~1.5 GB, not
in Git). Chromosomes are unprefixed (`1`, `4`, `X`). NCBI `COMMON=1`: at
least one 1000 Genomes population has MAF ≥ 1%, with ≥ 2 founders
contributing that allele. Frequency fields: `CAF` (first value = reference)
and `TOPMED`. There is no standard `AF=`.

```bash
gzip -dc data/dbSNP_b151_GRCh38p7_common_all_20180418.vcf.gz | head -n 20
```


| Header                  | Must be         |
| ----------------------- | --------------- |
| `##source=dbSNP`        | NCBI dbSNP      |
| `##dbSNP_BUILD_ID=151`  | build 151       |
| `##reference=GRCh38.p7` | patch 7         |
| `##fileDate=20180418`   | this extract    |
| `#CHROM` data           | `1`, not `chr1` |


Example FASTA headers use `chrom=4` and 0-based BED; the VCF uses `4` and
1-based `POS` (the reader subtracts 1). Do not mix GRCh37/hg19 or
later-patch / alt-contig coordinates.

## BLAST subject (GENCODE 45, GRCh38.p14, CHR)

[https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz)

Release notes: [https://www.gencodegenes.org/human/release_45.html](https://www.gencodegenes.org/human/release_45.html)

Local FASTA: `data/gencode.v45.transcripts.fa` (~454 MB). BLAST prefix:
`data/gencode_v45_transcripts_db`. Comprehensive CHR transcripts including
MT, **not** `pc_transcripts`. Queries are DNA (`T`, never `U`).

`bsst db init --gencode-v45-transcripts` builds the library. Equivalent
manual command (do **not** pass `-parse_seqids`; GENCODE headers contain
`|`):

```bash
makeblastdb \
  -in data/gencode.v45.transcripts.fa \
  -dbtype nucl \
  -out data/gencode_v45_transcripts_db \
  -title "GENCODE v45 transcripts CHR GRCh38.p14"
```

```bash
grep -c '^>' data/gencode.v45.transcripts.fa
head -n 1 data/gencode.v45.transcripts.fa
```


| Check                     | Requirement                                    |
| ------------------------- | ---------------------------------------------- |
| Record count              | **252930**                                     |
| First header              | `ENST00000456328.2` / `DDX11L2-202` / `lncRNA` |
| Last (MT) header          | `MT-TP-201` / `Mt_tRNA`                        |
| `protein_coding`          | 89110                                          |
| `nonsense_mediated_decay` | 21427                                          |
| `lncRNA`                  | 57722                                          |


GENCODE 47+ has ~385k transcripts and may still start with DDX11L2. Do not
identify the release from the first header alone.

BLAST compares **transcript sequence**, not genomic intervals, so p7 vs p14
does not shift VCF-style coordinates. The off-target subject is the
January 2024 GENCODE 45 set, not live Ensembl.