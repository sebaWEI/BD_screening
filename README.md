# Binding Site Safety Filter (`bssf`)

Dry-lab **binding-site safety filter** (BD_screening package) for
**iGEM PekingHSC 2026** (HEPHA-RNA).

It keeps sense sites on a transcript 3′UTR whose antisense element does
**not** have a same-strand BLAST hit in another transcript's 5′UTR or 3′UTR.
Optionally (`--variants`) it also drops sites that overlap high-frequency
common SNPs.

That is the whole job: offtarget UTR gate, optional population-variation
gate. A pass is not upregulation.

- [Team model repo](https://github.com/sebaWEI/PekingHSC-2026-Model)
- [iGEM team](https://teams.igem.org/6371)
- [Wiki](https://2026.igem.wiki/pekinghsc/)

Each run writes `runs/<timestamp>_<id>/`.

## Pipeline

1. 3′UTR from `--utr`, or Ensembl 111 canonical coding 3′UTR for `--gene`.
2. Sense sites from `--sites`, or 40 nt windows (step 1). A supplied site
   must occur exactly once in the UTR.
3. Drop when `blastn-short` finds a same-strand, non-self alignment ≥ 20 nt
   that overlaps another transcript's 5′UTR or 3′UTR. CDS-only hits stay.
   Raise the length with `--offtarget-min-length`. `--gene` is required so
   the intended gene is not called offtarget.
4. Optionally `--variants`: drop dbSNP COMMON (GRCh38.p14) overlaps.

Details: [docs/guide.md](docs/guide.md).

## Setup

Need Git, internet, and ~1 GB free for GENCODE. Variants use a ~3 MB
tabix index plus on-demand HTTPS queries (no 28 GB dbSNP download).

```bash
git clone https://github.com/sebaWEI/BD_screening.git
cd BD_screening
./scripts/setup.sh          # Windows: .\scripts\setup.ps1
```

Setup installs `uv` + `.venv`, BLAST+/tabix when available, downloads the
pinned DBs, builds the BLAST index and UTR table, then runs
`check_requirements`.

`check_requirements` must show `ok` for `blastn`, `makeblastdb`,
`blast_db`, `blast_db_identity`, `gencode_fasta`, `gencode_fasta_identity`,
and `gencode_utr_index`.

## Usage

```bash
source .venv/bin/activate   # Windows: .venv\Scripts\Activate.ps1
bssf filter --gene LETM1 --sites examples/LETM1.sites.fasta
bssf filter --gene LETM1 --sites examples/LETM1.sites.fasta --variants
```

Omit `--sites` and every 40 nt window goes into one `blastn` (300 s
timeout). Full LETM1/NSD2 UTRs will not finish in that limit—pass `--sites`.

| Output | Meaning |
|--------|---------|
| `candidates.tsv` | Sites with `status=pass` |
| `all_binding_sites.tsv` | Every site + `failure_reason` |
| `blast_matches.tsv` | Per site/gene offtarget summary; `drops_site` |
| `blast_hits.tsv` | Full BLAST table |
| `run.log` / `manifest.json` | Commands and versions |

Examples: `examples/LETM1.fasta`, `examples/NSD2.fasta`, and the matching
`*.sites.fasta` tiles.

| Command | Purpose |
|---------|---------|
| `bssf check_requirements` | PATH + pinned files |
| `bssf resources` | Producer names / URLs |
| `bssf db init --dbsnp-common-all --gencode-v45-transcripts` | Download DBs, BLAST, UTR index |
| `bssf filter --gene SYMBOL --sites FILE` | BLAST filter |
| `bssf filter --gene SYMBOL --sites FILE --variants` | BLAST + variants |
| `bssf config show` | Local config |

## If the check failed

Fix the `missing` / `mismatch` row, then re-run `bssf check_requirements`.
`variant_vcf` is only needed for `--variants`. `tabix` is optional.

```bash
# uv
curl -LsSf https://astral.sh/uv/install.sh | sh
uv python install 3.12 && uv sync

# BLAST+ / tabix
brew install blast htslib
# or: conda install -y -c bioconda blast htslib

# databases
bssf db init --dbsnp-common-all --gencode-v45-transcripts
```
