# Methods

Default: BLAST-filter sense sites on a transcript-oriented 3′UTR. Sites
come from `--sites`, or from 40 nt windows stepped by 1. `--gene` is
required so self hits are ignored. Variant overlap and RNAup ranking are
off unless requested. There is no composite score.

```text
3′UTR (file or Ensembl 111)
    → supplied sense sites, or windows (40 nt)
    → drop sites that are missing or repeated on the UTR
    → BLAST against GENCODE 45 transcripts
    → keep same-strand hits only; map them onto 5′UTR / 3′UTR
    → drop sites with a UTR hit ≥ --offtarget-min-length (default 20)
    → optional dbSNP overlap filter (--variants)
    → optional RNAup on the sites that passed (--rnaup)
    → rank those by ΔG_total
```

## Filters

**Variants.** NCBI dbSNP b151 `common_all_20180418` on GRCh38.p7
(chromosome names `1`, `4`, … not `chr1`). Frequencies are `CAF` / `TOPMED`,
not `AF=`. With no `--min-af`, every record in the extract is used. With
`--min-af`, the reader takes `AF` if present, else the largest
non-reference `CAF`; records with neither field are ignored. A tabix index
reads only the target 3′UTR; without it the same interval is kept, but the
~1.5 GB file is decompressed once.

`--variants` is off by default. The UTR needs `chrom`, 0-based BED
`start`/`end`, and strand. Incomplete coordinates skip variant filtering
for the whole run. Site coordinates are the UTR interval plus the site
offset.

**BLAST.** `blastn -task blastn-short` against GENCODE 45 CHR transcripts
(GRCh38.p14, 252930 sequences, all biotypes). Queries are DNA (`T`, never
`U`). A same-strand, non-self hit of **≥ 20 nt** that overlaps a 5′UTR or 3′UTR
flags the window. CDS-only hits and transcripts with no UTR annotation do
not. A hit with `sstart` > `send` is the element sequence on another
transcript, not a site the element can bind, and is ignored. The dropped
site's `failure_reason` names the gene, transcript, UTR, identity, and
both intervals. Identity and coverage
floors default to 0. Self hits match gene symbols and ENST accessions with
or without `.version`. Do not pass `-parse_seqids` to `makeblastdb`
(GENCODE headers contain `|`).

`blast_matches.tsv` is the table to read. It keeps the longest same-strand
alignment per site and gene. `region` is `5'UTR`, `3'UTR`, or both when
the alignment overlaps that part of the spliced transcript. `drops_site`
is `yes` only for those UTR hits that also meet the length cutoff.

All pending sites go out in one `blastn` process, which is stopped after
300 seconds. Use `--sites` for the sequences you intend to test. Omitting
`--sites` searches every 40 nt window of the UTR and does not finish for
LETM1 or NSD2 within that limit.

## RNAup

Each BD is scored against a target slice with 120 nt context on each side:

`--interaction_first --window N --temp T` (optional `--include_both`).

RNAup reports 1-based local coordinates; bsst maps them to 0-based UTR
offsets. The whole predicted target interaction must lie inside the site
(`anchor_overlap = 1`). Off-anchor hits are dropped. This section runs
only with `--rnaup`, and only on sites that already passed BLAST and,
when requested, the variant filter.

A long run of `(` / `)` in RNAup output is expected: the BD is the reverse
complement of the window, so the MFE is usually a ~40 bp intermolecular
helix. A `.` is an unpaired base (often a weak terminal A-U), not a crash.

## Run directory

Every invocation writes `runs/<UTC timestamp>_<run id>/`:

| File | Contents |
|------|----------|
| `manifest.json` | status, parameters, input SHA-256, tool versions, databases |
| `run.log` | stages plus full external commands, stdout, stderr |
| `inputs/target.fasta` | exact analyzed sequence |
| `all_binding_sites.tsv` | every site and `failure_reason` |
| `candidates.tsv` | sites with `status=pass`. With `--rnaup`, only anchored sites (`status=eligible`), ordered by ΔG_total |
| `blast_hits.tsv` | raw BLAST hits plus `is_self`, `match_region`, and `offtarget_risk` |
| `blast_matches.tsv` | longest same-strand hit per site and gene, with `region` and `drops_site` |

`status=pass` is the default result. Energy columns, `interaction_*`,
`anchor_overlap`, and `rank` are filled only when you passed `--rnaup`.
A dropped site's `failure_reason` looks like
`3'UTR of SIPA1L1 (SIPA1L1-207), 21 nt, 95.2% identity, site 20-40, transcript 1126-1146`.

| Field | Meaning |
|-------|---------|
| genomic `start`/`end` | 0-based half-open |
| `bd_sequence` | reverse complement of the sense site, DNA notation |
| `match_region` / `region` | `5'UTR`, `3'UTR`, `5'UTR+3'UTR`, or empty when the hit misses both |
| `drops_site` | `yes` when that UTR hit also meets `--offtarget-min-length` |
| `interaction_target_*`, `interaction_query_*` | RNAup 1-based inclusive, local (`--rnaup` only) |
| `interaction_utr_*` | mapped 0-based half-open UTR offsets (`--rnaup` only) |
| energy columns | kcal/mol (`--rnaup` only) |
| `anchor_overlap` | fraction of the predicted target interaction inside the site (`--rnaup` only) |
| `rank` | 1 = most negative ΔG_total (`--rnaup` only) |

If `candidates.tsv` is empty, read `all_binding_sites.tsv` → `failure_reason`:
`not_in_utr`, `ambiguous`, a 5′UTR or 3′UTR match description, and with `--variants`
`variant_overlap`. With `--rnaup`, a site that passed BLAST can still fail
as `off_anchor_interaction` or `RNAup_*`. Do not pool these tables with
runs from before the filter refactor.

Pinned files: [resources.md](resources.md).
What the scores are not: [limitations.md](limitations.md).

## References

Tool and database papers the pipeline uses. Titles are the publishers’
article titles.

- Mückstein et al. (2006). Thermodynamics of RNA–RNA binding. *Bioinformatics*
  **22**, 1177–1182. https://doi.org/10.1093/bioinformatics/btl024
- Lorenz et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular
  Biology* **6**, 26. https://doi.org/10.1186/1748-7188-6-26
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
