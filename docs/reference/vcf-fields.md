---
title: VCF fields
description: Every INFO and FORMAT tag GraffiTE writes, its meaning, and the code that sets it.
---

# VCF fields

!!! info "Applies to GraffiTE v1.1"
    This page documents the `v1.1dev` branch at commit `18a76d9`.
    Behaviour described in the [2024 paper](https://www.nature.com/articles/s41467-024-53294-2)
    corresponds to v1.0 and differs in places.


!!! warning "Partial"
    The HERV-K section below is written. The other sections are still outlines.


## Reading a GraffiTE VCF record


## Repeat annotation INFO fields


## TSD and polyA INFO fields


## HERV-K INFO fields

Written only under `--human`, and only to `pangenome.human.vcf` — the discovery
VCF `pangenome.vcf` is left untouched because it induces the graph. See
[HERV-K (HML-2) biology](../background/hervk-hml2.md) for what the states mean.

| field | |
|---|---|
| `HERVK_CLASS` | `null_solo`, `solo_prov`, `truncated_prov`, `null_prov`, `copy_number`, `other` |
| `HERVK_ALLELE_REF` | HERV-K state of the REF allele: `null`, `solo`, `provirus`, `prov_xN`, `partial`, `.` |
| `HERVK_ALLELE` | state of each ALT allele, in ALT order |
| `HERVK_EVIDENCE` | what resolved it: `ARCH_2LTR`, `ARCH_PERM`, `ARCH_INT_PERM`, `CNV_PERIOD`, `ARCH_SOLO`, `REF_ANNOT`, `DENOVO_LTR`, `UNRESOLVED`, `NON_HML2` |
| `HERVK_ARCH` | element architecture with consensus intervals, e.g. `LTR:575-968/INT:1-7536/LTR:1-574` |
| `HERVK_K` | LTR permutation point. An alignment property — see the caveat below |
| `HERVK_J` | internal-region permutation point, the `ARCH_INT_PERM` counterpart of `HERVK_K` |
| `HERVK_N_UNITS_REF` | proviral units in the reference element, a junction LTR counted once |
| `HERVK_UNIT_BP` | period of the reference array: one internal region plus one LTR |
| `HERVK_REF_STATE` | state of the masked reference window |
| `HERVK_LAMBDA` / `HERVK_NU` | bp of HML-2 LTR / internal sequence on the variant allele |
| `HERVK_COV` | fraction of the variant that is HML-2 sequence |
| `HERVK_PMAP` | confidence in the resolved class. **Reporting only** — it does not decide the class |
| `HERVK_NOTE` | `REF_ARCH_CONFLICT`, `CNV_UNITS:a->b`, `UNIT_COUNT_ASSUMED`, `INS_INTO_NONEMPTY_REF` |
| `HERVK_LOCUS`, `HERVK_LOCUS_N`, `HERVK_MERGE_FLAG` | locus grouping; `MERGE_FLAG` marks a locus holding more than one record |

After consolidation (`4_Genotyping/GraffiTE.merged.genotypes.human.vcf.gz`):

| field | |
|---|---|
| `HERVK_AC` / `HERVK_AN` | allele count per ALT, and total alleles called |
| `HERVK_AC_DISC` / `HERVK_AN_DISC` | the same from the assembly callset, as an independent check |
| `HERVK_DISC_CONCORDANT` | set when the two agree on every count |
| `HERVK_N_RESOLVED`, `HERVK_N_PARTIAL`, `HERVK_N_PLOIDY_EXCEEDED` | how the genotypes resolved |
| `HERVK_MEMBERS`, `HERVK_MEMBERS_MASKED` | records the locus consolidates |

!!! warning "Three things to read carefully"
    **Quote `HERVK_AC` and `HERVK_AN`, not `AF` alone.** `AN` below `2N` means
    members were structurally uncalled for some samples — the graph can find
    every carrier and still fail to confirm the non-carriers, which deflates
    `AN` without touching `AC`.

    **`HERVK_N_PLOIDY_EXCEEDED > 0`** means a sample's summed allele dosage
    exceeded its ploidy: a third allele lost to `bcftools norm -m-`. Those
    genotypes are set missing rather than guessed at.

    **`HERVK_K` is not a biological measurement.** It is where the aligner
    happened to split the reference LTR, it varies between haplotypes and
    callers, and nothing should key on its value.


## FORMAT fields


## Fields inherited from upstream callers


## Retired fields


## Filtering GraffiTE VCFs with bcftools
