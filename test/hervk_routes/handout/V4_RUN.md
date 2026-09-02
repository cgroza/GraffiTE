# v4 run: consolidated human VCFs

Read `README_HPC_CLAUDE.md` for the setup and `E2E_TEST.md` for stage E. V3_RUN.md
covers the copy-number work, which is merged. This file covers only what is new.

Pipeline is at **`feat/hervk-consolidated-vcfs`**, branched off `v1.1dev` after
PR #96 landed. Seven commits. `INPUTS.env` already names it, so `./bootstrap.sh`
pulls the right thing.

## What changed

A user got a HERV-K locus as a scatter of records, with nothing in either
delivered VCF to say which records belonged together or whether the locus was
an insertion polymorphism. Both VCFs now carry a merged view.

| # | change | changes output? |
|---|---|---|
| 1 | `build_locus_record` splices the REF span from the FASTA when no member spans the locus | **yes**, chr7-shaped loci consolidate |
| 2 | graph genotypes are withheld per allele, not per locus | **yes**, chr6's solo allele keeps its counts |
| 3 | `HERVK_MEI`, `HERVK_SOLO_PROV`, `HERVK_CNV` replace one exclusive category | **yes**, new fields |
| 4 | `pangenome.human.consolidated.vcf` | new file |
| 5 | a locus the `--human` filter cut in half is annotated, not dropped | **yes**, chr7 and chr8 keep their locus id |
| 6 | `hervk_reconcile` transfers the HERV-K INFO onto the genotyped subset | **yes**, 31 loci annotated instead of 4 |

Change 6 is the one to look at first. `merge_VCFs` transfers INFO from
`pangenome.vcf`, which carries no HERV-K annotation, so in the last run only 4
of 29 HERV-K records in `GraffiTE.merged.genotypes.human.vcf.gz` had a locus
id, and you could not filter any of them on `HERVK_MEI`.

## Cost

Same as the v3 run. Nothing here re-masks or re-genotypes. The two consolidation
steps read VCFs and make one `samtools faidx` call per multi-allelic locus.

## Run it

Unchanged:

```bash
./bootstrap.sh          # pulls feat/hervk-consolidated-vcfs
$EDITOR INPUTS.env      # paths, as before
./preflight.sh
./run_hervk_test.sh
```

Then, from the handout directory:

```bash
python3 assert_hervk_test.py --outdir <your --out directory>
```

## Nextflow was never executed against these changes

Nextflow is not installed on the machine these commits were written on, so
nobody has parse-checked the Groovy. I rendered both process script blocks to
plain bash and they pass `bash -n`, and I ran the commands they produce by hand
against the v3 outputs. Your launch is the first real test of the wiring.

Two places to check if a process dies early:

- `hervk_annotate` calls `hervk_reconcile.py consolidate --source discovery`
  with `--reference "$REF"`, the same `$REF` the masking steps above it use.
- `hervk_reconcile` bgzips and indexes `merged.human.vcf` before
  `bcftools annotate` reads it. htslib needs the index even though the file is
  only streamed.

## What must be true

`assert_hervk_test.py` checks every row below and exits non-zero on any of
them. The counts come from re-running the current code over the v3 outputs, not
from a fresh pipeline run, so treat a small difference as something to look at
rather than a failure.

| check | expected |
|---|---|
| `3_TSD_search/pangenome.human.consolidated.vcf` | exists |
| its record count | below `pangenome.human.vcf`; that file keeps every record, because it induces the graph |
| HERV-K loci in each consolidated VCF | 31 |
| `HERVK_MEI` | 21 |
| `HERVK_SOLO_PROV` | 9 |
| `HERVK_CNV` | 3 |
| `HERVK_LOCUS_INCOMPLETE` | 3 |
| the two consolidated VCFs | identical on all five counts |
| chr6:78,894,316 | carries `HERVK_SOLO_PROV` and `HERVK_CNV`, and not `HERVK_MEI` |
| chr6 `HERVK_ALLELE_NOGT` | `prov_x2` |
| chr6 `HERVK_AC` | starts `8,` |
| chr6 `HERVK_AC_DISC` | `8,1` |
| chr7:4,699,714 and chr8:7,552,031 | carry a locus id, `HERVK_LOCUS_INCOMPLETE` and `HERVK_ALLELE_SET` |
| `hervk_loci.tsv` | has `locus_type`, `mei`, `solo_prov`, `cnv` columns |
| every flag | declared in both VCF headers |
| chr11:101,704,640 | still AC=23, AN=40, `HERVK_DISC_CONCORDANT` |

Do not skip the header check. An undeclared flag makes bcftools assume
`Type=String`, and `-i 'INFO/HERVK_MEI=1'` then matches nothing and reports no
error.

## What will move, and should

`GraffiTE.merged.genotypes.human.vcf.gz` loses records to consolidation.
`hervk_unconsolidated_records.vcf` gains the members that were merged, plus our
INFO definitions in its header.

`HERVK_AN` at chr6 is 22, not 40. `vg call` emits no call for 9 of the 20
samples at that record, and consolidation only reports what it found.
`HERVK_AN_DISC=40` sits beside it.

## Report back

`bundle_results.sh` as before, plus:

- both consolidated VCFs
- `3_TSD_search/hervk_discovery_consolidation_report.md`
- the full `assert_hervk_test.py` output, pass or fail
