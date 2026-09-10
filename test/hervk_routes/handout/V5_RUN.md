# v5 run: the chr7 locus, whole

Read `README_HPC_CLAUDE.md` for the setup. V4_RUN.md covers the consolidated
VCFs, which are merged. This file covers one parameter change.

Pipeline is at **`fix/hervk-pair-rule-n-hits`**, branched off `v1.1dev` after
PR #97. One commit.

## What changed

`--hervk_pair_max_hits`, new, default 3. The `--human` filter's HERV-K carve-out
tested `n_hits==2`; it now tests `n_hits<=3`.

`SVA_A` shares homology with `LTR5_Hs`, so RepeatMasker splits a small SVA hit
off the LTR of an HML-2 provirus, which is what the carve-out is for. It splits
the internal region of a degraded or rearranged provirus too, giving a third
hit. At chr7:4,699,714 two of the three records come back as
`LTR5_Hs,SVA_A,HERVK-int` and `HERVK-int,SVA_A,HERVK-int`, and the equality test
dropped both. The locus kept one record and reported 39/1 for a two-unit against
three-unit difference where the truth is 26/12/2 across three states.

Nothing else in the clause moved. It still requires `LTR/ERVK` beside
`Retroposon/SVA`, an `^HERVK-int` id, and `|SVLEN|` no greater than one
provirus, and those three are what keep the exception narrow.

## Measured before implementing

Against `3_TSD_search/pangenome.vcf` from the June run, 132,511 records.
Applying the shipped filter to it with bcftools returns 5,817, matching the
delivered `pangenome.human.vcf`, so the two runs share an input and the deltas
below are comparable.

| filter | kept | Alu | L1 | SVA | VNTR | HERV-K |
|---|---|---|---|---|---|---|
| v1.1 as shipped | 5,817 | 4,874 | 701 | 114 | 93 | 35 |
| **pair rule `n_hits<=3`** | **5,819** | 4,874 | 701 | 114 | 93 | **37** |
| `n_hits>=1` for LTR/ERVK | 5,859 | 4,874 | 701 | 114 | 93 | 77 |
| `n_hits<=3` in the single-hit rule | 5,837 | 4,874 | 701 | 114 | 93 | 55 |

No non-HERV-K count moves under any of them. The pair-rule change admits
exactly the two chr7 records and nothing else. The two blunt variants admit 18
and 40 more that are an LTR5 fragment beside something unrelated:
`COMP-subunit_FAM90A,LTR5A` eight times, `ALR/Alpha,LTR5A,ALR/Alpha`,
`L1PA10,LTR5_Hs`, `FLAM_C,MER67B,LTR5_Hs`, and `HERVK9-int(x),AluYm1`, a
lineage `human_hervk_ids` deliberately excludes.

## What this does not fix

Two loci still lose a member, and neither loses an allele state:

| locus | absent member | blocked by |
|---|---|---|
| chr1:75,219,429 | `chr1-75220275-INS-16215` | the pair rule's `abs(SVLEN)<=10500`; it is 16,215 bp |
| chr8:7,552,031 | `chr8-7552072-INS-218` | the global `abs(SVLEN)>=250` floor; it is 218 bp |

Both keep `HERVK_LOCUS_INCOMPLETE`. chr8's would need the size floor lowered for
every class, which is a different question.

## Run it

Unchanged from v4, except the branch:

```bash
REVISION=fix/hervk-pair-rule-n-hits ./bootstrap.sh
$EDITOR INPUTS.env      # set REVISION and OUTDIR; keep the v4 paths
./preflight.sh
./run_hervk_test.sh
```

Set `GENOTYPED_VCF` as before. Stage E is half of what this tests.

## What must be true

`assert_hervk_test.py` checks all of it and exits non-zero on any of it.

| check | expected |
|---|---|
| `pangenome.human.vcf` record count | **5,819**, up from 5,817 |
| chr7:4,699,714 | one locus of 3, no `HERVK_LOCUS_INCOMPLETE` |
| chr7 alleles | `prov_x2`, `prov_x3` and `provirus` |
| chr7 allele counts | 12, 2 and 26 of 40 |
| `HERVK_LOCUS_INCOMPLETE` | **2**, on chr1:75,219,429 and chr8:7,552,031 |
| HERV-K loci | 31, unchanged |
| `HERVK_MEI` / `HERVK_SOLO_PROV` / `HERVK_CNV` | 21 / 9 / 3, unchanged |
| Alu, L1, SVA, SVA VNTR-only | 4,874 / 701 / 114 / 93 at stage 1, unchanged |
| consolidated loci at stage E | 5, up from 4 |

chr7 has three members in the subset now, so the reconciler merges them into
one multi-allelic record instead of annotating a lone survivor. All three are
copy-number records whose graph genotypes the pipeline withholds, so the locus
reports `HERVK_ALLELE_NOGT` for its alleles and takes its counts from the
discovery genotypes. chr12:133,148,144 already exercises that path.

## Simulated before asking for the cluster

I spliced the two records into the delivered v4 `pangenome.human.vcf` and ran
`hervk_reconcile.py flag` and `consolidate` over the result with the code on this
branch. chr7 consolidated to one record reporting 26/12/2, the split-locus count
fell from 3 to 2, and nothing else in the locus layer moved. The Nextflow wiring
is the part that has not run, as usual.

## Report back

`bundle_results.sh`, plus both consolidated VCFs, the new
`human_filter_summary.txt`, and the full assertion output.
