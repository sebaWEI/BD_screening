# Methods

Default: BLAST-filter sense sites on a transcript-oriented 3′UTR. Sites
come from `--sites`, or from 40 nt windows stepped by 1. `--gene` is
required so self hits are ignored. A site is dropped only for a same-strand
hit of at least `--offtarget-min-length` nt (default 20) that overlaps a
5′UTR or 3′UTR. Variant overlap (`--variants`) is off unless requested.
There is no energy score and no composite rank.

```text
3′UTR (file or Ensembl 111)
    → supplied sense sites, or windows (40 nt)
    → drop sites that are missing or repeated on the UTR
    → BLAST against GENCODE 45 transcripts
    → keep same-strand hits only; map them onto 5′UTR / 3′UTR
    → drop sites with a UTR hit ≥ --offtarget-min-length (default 20)
    → optional dbSNP overlap filter (--variants)
    → write pass/fail tables
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
`U`). A same-strand, non-self hit of at least `--offtarget-min-length` nt
(default **20**) that overlaps a 5′UTR or 3′UTR flags the window. CDS-only
hits and transcripts with no UTR annotation do not. A hit with `sstart` >
`send` is the element sequence on another transcript, not a site the
element can bind, and is ignored. The dropped site's `failure_reason`
names the gene, transcript, UTR, identity, and both intervals. Identity
and coverage floors default to 0. Self hits match gene symbols and ENST
accessions with or without `.version`. Do not pass `-parse_seqids` to
`makeblastdb` (GENCODE headers contain `|`).

`blast_matches.tsv` is the table to read. It keeps one same-strand
alignment per site and gene: a UTR hit that drops the site when present,
otherwise the longest same-strand hit. `region` is `5'UTR`, `3'UTR`, or
both when the alignment overlaps that part of the spliced transcript.
`drops_site` is `yes` only for those UTR hits that also meet the length
cutoff.

All pending sites go out in one `blastn` process, which is stopped after
300 seconds. Use `--sites` for the sequences you intend to test. Omitting
`--sites` searches every 40 nt window of the UTR and does not finish for
LETM1 or NSD2 within that limit.

## Run directory

Every invocation writes `runs/<UTC timestamp>_<run id>/`:

| File | Contents |
|------|----------|
| `manifest.json` | status, parameters, input SHA-256, tool versions, databases |
| `run.log` | stages plus full external commands, stdout, stderr |
| `inputs/target.fasta` | exact analyzed sequence |
| `all_binding_sites.tsv` | every site and `failure_reason` |
| `candidates.tsv` | sites with `status=pass` |
| `blast_hits.tsv` | raw BLAST hits plus `is_self`, `match_region`, and `offtarget_risk` |
| `blast_matches.tsv` | longest same-strand hit per site and gene, with `region` and `drops_site` |

`status=pass` is the default success label. A dropped site's
`failure_reason` looks like
`3'UTR of SIPA1L1 (SIPA1L1-207), 21 nt, 95.2% identity, site 20-40, transcript 1126-1146`.

| Field | Meaning |
|-------|---------|
| genomic `start`/`end` | 0-based half-open |
| `bd_sequence` | reverse complement of the sense site, DNA notation |
| `match_region` / `region` | `5'UTR`, `3'UTR`, `5'UTR+3'UTR`, or empty when the hit misses both |
| `drops_site` | `yes` when that UTR hit also meets `--offtarget-min-length` |

If `candidates.tsv` is empty, read `all_binding_sites.tsv` → `failure_reason`:
`not_in_utr`, `ambiguous`, a 5′UTR or 3′UTR match description, and with `--variants`
`variant_overlap`.

Pinned files: [resources.md](resources.md).
What a pass is not: [limitations.md](limitations.md).

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
