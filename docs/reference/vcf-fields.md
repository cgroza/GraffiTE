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

Written only under `--human`. `pangenome.vcf` is left untouched because it
induces the graph, so these fields appear in `pangenome.human.vcf`, in the two
consolidated VCFs described below, and in `hervk_candidates.vcf`. See
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

### What varies at the locus

Three independent flags, any of which may be set. A locus can be more than one
of these at once, which is why they are separate fields.

| field | |
|---|---|
| `HERVK_MEI` | a `null` allele segregates: the element is absent from some haplotypes, so a transposition produced the difference. **Filter on this one.** It separates an insertion polymorphism, comparable to an *Alu*, L1 or SVA insertion, from structural variation in an element every haplotype carries |
| `HERVK_SOLO_PROV` | a solo LTR and a provirus both segregate |
| `HERVK_CNV` | some allele carries two or more proviral units |
| `HERVK_LOCUS_TYPE` | a single summary label: `null_vs_present`, `copy_number`, `solo_vs_provirus`, `unresolved` |

GraffiTE derives `HERVK_LOCUS_TYPE` from the flags, so one label has to stand
in for all three and a locus that is two things loses half of what it is. In
the CaG cohort chr6:78,894,316 segregates a solo LTR, a provirus and a two-unit
allele. It is labelled `copy_number` and carries both `HERVK_SOLO_PROV` and
`HERVK_CNV`. Nine loci carry `HERVK_SOLO_PROV`; seven are labelled
`solo_vs_provirus`.

### Consolidated records

GraffiTE can emit several VCF records for one HERV-K locus: two breakpoints
for the same insertion, or a deletion and an insertion describing opposite
directions of the same event. Consolidating collapses them onto one
multi-allelic record.

This happens at both stages, and each time GraffiTE writes a new file rather
than changing the old one. It cannot change the inputs. The discovery VCF
induces the graph, and the graph VCF is the native record of what `vg call`
did.

| this file | is | holds |
|---|---|---|
| `3_TSD_search/pangenome.human.vcf` | unmerged | one record per call, as the graph needs it |
| `3_TSD_search/pangenome.human.consolidated.vcf` | merged | one record per locus, genotypes from the assemblies |
| `4_Genotyping/GraffiTE.merged.genotypes.vcf.gz` | unmerged | every genotyped call, `--human` or not |
| `4_Genotyping/GraffiTE.merged.genotypes.human.vcf.gz` | merged | one record per locus, genotypes from the graph |
| `4_Genotyping/hervk_unconsolidated_records.vcf` | the members | what each merged record was built from |

The `##hervk_consolidation` header line names the source.

| field | |
|---|---|
| `HERVK_AC` / `HERVK_AN` | allele count per ALT, and total alleles called |
| `HERVK_AC_DISC` / `HERVK_AN_DISC` | the same from the assembly callset, as an independent check |
| `HERVK_DISC_CONCORDANT` | set when the two agree on every count |
| `HERVK_DISC_PLOIDY_MISMATCH` | the two callsets disagree on ploidy here, so GraffiTE withholds the assembly counts rather than report them on a denominator the two do not share. Expect it on hemizygous chromosomes and with a haploid discovery caller such as SVIM-asm run per haplotype |
| `HERVK_N_RESOLVED`, `HERVK_N_PARTIAL`, `HERVK_N_PLOIDY_EXCEEDED` | how the genotypes resolved |
| `HERVK_MEMBERS`, `HERVK_MEMBERS_MASKED` | records the locus consolidates |
| `HERVK_ALLELE_NOGT` | alleles `HERVK_AC` reports as 0 because the graph cannot count them, not because they are absent. Their counts are in `HERVK_AC_DISC` |
| `HERVK_ALLELE_SET` | every allele state at the locus, including states carried by records that are not in this file |
| `HERVK_MEMBERS_ABSENT`, `HERVK_LOCUS_INCOMPLETE` | the `--human` filter removed a member, so the counts here do not cover every allele |

!!! warning "A copy-number allele cannot be genotyped from short reads"

    Its ALT path repeats sequence the REF path already carries, so a read
    from the pre-existing copy traverses it and non-carriers pick up ALT
    support. More depth does not remove that. GraffiTE withholds those alleles
    from the graph genotypes and names them in `HERVK_ALLELE_NOGT`. The other
    alleles at the same locus are genotyped normally, and the withheld one
    takes its frequency from the assemblies. At chr6:78,894,316 the graph
    counts the solo allele at AC=8, matching the assemblies, and the two-unit
    allele comes from `HERVK_AC_DISC`.

!!! warning "`HERVK_LOCUS_INCOMPLETE` means the frequencies are partial"

    The `--human` filter is narrower than the HERV-K candidate list, so it can
    remove some members of a locus and keep others. Three CaG loci are in that
    state. At chr7:4,699,714 it keeps one of three. All three describe the
    same ~8.5 kb unit, and the two it removes carry a third RepeatMasker
    fragment that the filter's HERV-K clause does not admit. The surviving
    record
    carries the locus id and the full `HERVK_ALLELE_SET`, but its counts cover
    only the alleles still present. `hervk_candidates.vcf` has all the records
    with their assembly genotypes.

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

Read HERV-K from a consolidated VCF. In the unmerged ones a locus is several
records and the same allele can be counted twice.

Insertion polymorphisms only, which is what to use alongside *Alu*, L1 and SVA:

```bash
bcftools view -i 'INFO/HERVK_MEI=1' GraffiTE.merged.genotypes.human.vcf.gz
```

Everything else HERV-K, the loci where the element is in every haplotype and
its structure varies:

```bash
bcftools view -i 'INFO/HERVK_LOCUS!="." && INFO/HERVK_MEI=0' \
    GraffiTE.merged.genotypes.human.vcf.gz
```

Drop loci whose counts do not cover every allele, either because the graph
could not count one or because the `--human` filter removed a member:

```bash
bcftools view -e 'INFO/HERVK_ALLELE_NOGT!="." || INFO/HERVK_LOCUS_INCOMPLETE=1' \
    GraffiTE.merged.genotypes.human.vcf.gz
```

SVA VNTR-only records, which RepeatMasker annotates as `Simple_repeat` rather
than `Retroposon/SVA`, are excluded with:

```bash
bcftools view -e 'INFO/repeat_ids~"VNTR_only"' pangenome.human.vcf
```
