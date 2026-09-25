# Binding Site Selection Tool

`bsst` (Binding Site Selection Tool) is the dry-lab **binding-site** filter
for **iGEM PekingHSC 2026** (HEPHA-RNA). It keeps sense binding sites on a
transcript 3′UTR whose antisense element does not have a same-strand BLAST hit
in another transcript's 5′UTR or 3′UTR. Optionally, you can use it
to filter out binding sites that overlap high-frequency variants. It also provides
an energy-calculation interface that gives some insight into those binding sites.

- [Main model repository of our team](https://github.com/sebaWEI/PekingHSC-2026-Model)
- [Team iGEM home page](https://teams.igem.org/6371)
- [Wiki home page](https://2026.igem.wiki/pekinghsc/)
- [Wiki model page](https://2026.igem.wiki/pekinghsc/model)
- [Wiki tutorial page](https://2026.igem.wiki/pekinghsc/documents)

Each run writes `runs/<timestamp>_<id>/` with the input, parameters, logs,
and pass/fail tables.

Hepha is specified to bind the **3′UTR**. A BLAST pass is not evidence of
up-regulation: [docs/limitations.md](docs/limitations.md).

## What a run does

1. Take a 3′UTR from `--utr`, or fetch the canonical coding 3′UTR for
   `--gene` from Ensembl 111.
2. Take sense sites from `--sites`. If `--sites` is not provided, slide 40 nt
   windows along that UTR. A supplied site must occur exactly once in the UTR.
3. Drop a site when `blastn-short` finds a same-strand, non-self alignment of
   at least 20 nt that overlaps another transcript's 5′UTR or 3′UTR. CDS-only
   hits are kept. Raise the length with `--offtarget-min-length`. `--gene` is
   required so the intended gene is not called an off-target.
4. With `--variants`, also drop dbSNP `common_all` overlaps. With `--rnaup`,
   score only the sites that passed and rank them by total ΔG.

Methods and output fields: [docs/methods.md](docs/methods.md).
Pinned files: [docs/resources.md](docs/resources.md).

## Setup

**Need:** Git, internet (GitHub, NCBI, EBI), about **3 GB** free disk, and
20–40 minutes the first time (most of that is a ~2 GB download).

### 1. Install Git and clone

- **macOS:** open **Terminal**. If `git --version` fails: `xcode-select --install`
- **Windows:** install [Git for Windows](https://git-scm.com/download/win), then open **PowerShell**.
- **Linux (Debian/Ubuntu):** `sudo apt update && sudo apt install -y git curl`

```bash
git clone https://github.com/sebaWEI/binding_site_selection_tool.git
cd binding_site_selection_tool
```

SSH alternative: `git clone git@github.com:sebaWEI/binding_site_selection_tool.git`

### 2. Run the setup script

The script installs `uv` and the Python project into `.venv`, installs
RNAup / BLAST+ / tabix with Homebrew when `brew` is on `PATH`, downloads
the pinned databases and GENCODE annotation, builds the BLAST index and the
spliced UTR table, and runs `check_requirements`.

**macOS / Linux:**

```bash
./scripts/setup.sh
```

**Windows (PowerShell):**

```powershell
.\scripts\setup.ps1
```

With Homebrew on `PATH`, this path does not use conda. Later `bsst` commands
use `.venv`, as in [Usage](#usage). The script does not install Git or Homebrew.
Machines without Homebrew, including Windows, can install the same native tools
with conda; that is only a fallback, in [If the check failed](#if-the-check-failed).

### 3. If the check passed, use bsst

`check_requirements` exits 0 when `blastn`, `makeblastdb`, `blast_db`,
`blast_db_identity`, `gencode_fasta`, `gencode_fasta_identity`, and
`gencode_utr_index` are `ok`.
Go to [Usage](#usage).

If any row is `missing` or `mismatch`, follow
[If the check failed](#if-the-check-failed).

## Usage

From the repo root, activate the project environment once. After that, type
`bsst` directly and leave off `uv run`:

```bash
source .venv/bin/activate
```

Windows PowerShell: `.venv\Scripts\Activate.ps1`.

The usual command names the gene and the sense sites. Each site must occur
exactly once in the UTR. Omitting `--utr` fetches the canonical coding 3′UTR
from Ensembl REST archive **111** (`https://e111.rest.ensembl.org`).

```bash
uv run bsst filter --gene LETM1 --sites examples/LETM1.sites.fasta
```

With the environment already active, that is `bsst filter --gene LETM1 --sites examples/LETM1.sites.fasta`.

Omit `--sites` and every 40 nt window (step 1) is sent in **one** `blastn`
process. That process stops after 300 seconds. A full LETM1 or NSD2 3′UTR
is about 2,900–3,200 windows and does not finish in that limit, which is
what `bsst filter --gene LETM1` hits.

That fetch never falls back to live Ensembl. Pass `--utr FILE` to use a
local 3′UTR instead, for example `examples/LETM1.fasta`. Fetching still
slides the whole UTR unless you also pass `--sites`.

`--variants` adds the dbSNP overlap filter. `--rnaup` scores only sites
that passed and ranks them by total ΔG (more negative first).
`--offtarget-min-length N` changes the UTR-hit length gate (default 20).

When it finishes it prints `Run completed: runs/<timestamp>_<id>/`.
Open that folder:

| File | What to look at |
| -------------------- | ------------------------------------------- |
| `candidates.tsv` | Sites that passed. `rank` and the energy columns stay empty unless you passed `--rnaup` |
| `all_binding_sites.tsv` | Every site. A dropped site's `failure_reason` names the gene, transcript, 5′UTR or 3′UTR, identity, and both intervals |
| `blast_matches.tsv` | One row per site and off-target gene. `region` is the UTR hit; `drops_site` is `yes` only when that hit also meets the length cutoff |
| `blast_hits.tsv` | Full BLAST table, including reverse-strand hits that are ignored for filtering |
| `run.log` | Commands that were actually executed |
| `manifest.json` | Tool versions and which databases were used |

The bundled examples are `examples/LETM1.fasta` / `examples/NSD2.FASTA`
and the wet-lab sense sites in `examples/LETM1.sites.fasta` (9 tiles) and
`examples/NSD2.sites.fasta` (18 tiles). Those site files are the reverse
complements of the SnapGene binding domains, so each sequence matches the
corresponding 3′UTR exactly once.

Commands below use `uv run` so they work before activation. After
`source .venv/bin/activate`, drop that prefix.


| Command                                                            | When you type it                          |
| ------------------------------------------------------------------ | ----------------------------------------- |
| `uv run bsst --help`                                               | Confirm the CLI is installed              |
| `uv run bsst check_requirements`                                   | PATH + pinned files after setup           |
| `uv run bsst resources`                                            | Print producer names and download URLs    |
| `uv run bsst db init --dbsnp-common-all --gencode-v45-transcripts` | Download the databases, build BLAST, and project 5′UTR / 3′UTR coordinates |
| `uv run bsst filter --gene SYMBOL --sites FILE`                    | Usual run: gene plus sense sites          |
| `uv run bsst filter --gene SYMBOL --sites FILE --offtarget-min-length 25` | Same run with a higher UTR-hit length gate |
| `uv run bsst filter --gene SYMBOL --utr FILE`                      | Slide the whole UTR; often hits the 300 s BLAST limit |
| `uv run bsst filter --gene SYMBOL`                                 | Fetch the 3′UTR, then slide it            |
| `uv run bsst config show`                                          | Show local config paths                   |

## If the check failed

Match the `missing` or `mismatch` row, run that command, then
`uv run bsst check_requirements` again. URLs: `uv run bsst resources`.
`RNAup` matters only for `--rnaup`. `variant_vcf` matters only for `--variants`.
`tabix` is optional.

**`uv`**

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
uv python install 3.12
uv sync
```

```powershell
# Windows
irm https://astral.sh/uv/install.ps1 | iex
uv python install 3.12
uv sync
```

**`RNAup`, `blastn`, `makeblastdb`, `tabix`**

```bash
# macOS / Linux with Homebrew
brew tap brewsci/bio
brew install brewsci/bio/viennarna blast htslib
```

```bash
# only when Homebrew is not available; this does not replace .venv
conda install -y -c conda-forge -c bioconda viennarna blast htslib
```

ViennaRNA builds: https://www.tbi.univie.ac.at/RNA/#download

**`blast_db`, `gencode_fasta`, `gencode_utr_index`, `variant_vcf`**

Do not download these in a browser. This command checks checksums, writes the
dbSNP VCF, GENCODE FASTA, and GENCODE GFF3 under `data/`, builds the BLAST
database, and writes `data/gencode.v45.utr_on_transcript.tsv`. That table is
how a transcript-coordinate BLAST hit is called 5′UTR or 3′UTR. `filter`
will not run without it:

```bash
uv run bsst db init --dbsnp-common-all --gencode-v45-transcripts
```

A checksum `mismatch` means the file is wrong or truncated. For the
GENCODE GFF3, `db init` deletes a bad copy and downloads again. Use
`--force` only when you intend to replace the local copies. If the FASTA is
already present but `blast_db` is not, build it yourself. Do **not** add
`-parse_seqids` (GENCODE headers contain `|`):

```bash
makeblastdb \
  -in data/gencode.v45.transcripts.fa \
  -dbtype nucl \
  -out data/gencode_v45_transcripts_db \
  -title "GENCODE v45 transcripts CHR GRCh38.p14"
```

Further reading: [docs/methods.md](docs/methods.md) ·
[docs/resources.md](docs/resources.md) ·
[docs/limitations.md](docs/limitations.md).

