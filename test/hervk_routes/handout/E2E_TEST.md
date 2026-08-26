# End-to-end test without re-genotyping

Graph genotyping 20 samples is the expensive stage, and nothing in the HERV-K
work depends on *how* the calls were made — only on having them. So stage E can
be pointed at a genotyped VCF from an earlier run, and the whole wiring runs
end to end for the cost of the discovery stage alone.

Set two things in `INPUTS.env`:

```bash
RM_DIR="/abs/path/2_Repeat_Filtering"                              # skips RepeatMasker
GENOTYPED_VCF="/abs/path/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz"
```

then `./preflight.sh && ./run_hervk_test.sh` as usual. The run adds
`--hervk_reconcile_vcf`, keeps `--genotype false`, and exercises:

- `hervk_annotate` — architecture, reference masking, classification, locus flags
- the `3_TSD_search` → `4_Genotyping` channel wiring
- `hervk_reconcile` — ID subsetting, dosage resolution, consolidation

Everything except `vg call` itself. Expected cost is re-run 2's ~34 minutes plus
a couple of minutes for stage E.

Sex correction does not matter here. The reconciler reads ploidy from the GT
field it is given and never imposes one, so an uncorrected VCF consolidates the
same way — it just carries whatever ploidy `vg call` wrote.

## New outputs in `4_Genotyping/`

```
GraffiTE.merged.genotypes.human.vcf.gz    human subset, HERV-K loci consolidated
hervk_unconsolidated_records.vcf          the member records removed, verbatim
hervk_reconciliation_report.md            per-locus AC/AN, and why
```

## What to check

Measured locally against the CaG data, so these are expectations, not guesses:

| check | expected |
|---|---|
| loci consolidated | 3 (chr1, chr11, chr12) |
| loci annotated in place | 1 (chr6 — its partner is a masked tandem) |
| loci skipped | 1 (chr8 — one member's ID did not survive `merge_VCFs`) |
| chr11 | `AC=23 AN=40`, and `HERVK_DISC_CONCORDANT` set |
| chr12 | `AC=6,21 AN=33`, `HERVK_N_PLOIDY_EXCEEDED=2` |
| chr6 tandem record | `HERVK_GT_MASKED`, every genotype `./.` |
| record accounting | in − 6 archived + 3 consolidated = out |

The chr12 ploidy violations are the point, not a fault: two samples whose summed
allele dosage exceeds ploidy because `bcftools norm -m-` flattened the third
allele. They are set missing rather than guessed at, and they are exactly the
gap between the graph and the assemblies at that locus.

`AN` below `2N` is likewise expected at chr6 and chr12 — members were
structurally uncalled (no depth reported at all) for some samples. That
depresses AN without touching AC: at chr6 the graph recovers every solo carrier
and only fails to confirm non-carriers. Report both numbers; `AF` alone hides it.

## If the genotyped VCF is from an older candidate set

Loci whose members were filtered out before that genotyping run will be skipped
with `N of M members present in the genotyped VCF`. That is correct behaviour —
the reconciler refuses to emit a partial locus — but it means the numbers above
only hold when the genotyped VCF and the discovery run share a candidate set.
chr8 is exactly this case in the CaG data.
