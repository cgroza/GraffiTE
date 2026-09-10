---
title: VCF fields
description: Every INFO and FORMAT tag GraffiTE writes, its meaning, and the code that sets it.
---

# VCF fields

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

Every table below gives the `Number` and `Type` from the header line the code writes, and the
**Source** column names the line that writes it. Fields are grouped by the stage that adds them.
Which files carry which group is in [Output files](outputs.md).

## Reading a GraffiTE VCF record

GraffiTE writes standard VCF 4.2 with sequence-resolved alleles. One record from the
`pangenome.human.vcf` of a three-sample PAV run, the ALT sequence and the PAV-specific INFO
fields shortened:

```text
chr1  18081  chr1-18082-INS-315_10  t  tAGAAGGAATAAGACGGGCCGGGT...  .  PASS
      ID=chr1-18082-INS-315;SVTYPE=INS;SVLEN=315;HAP=h2;...;NumCollapsed=4;CollapseId=5.0;
      n_hits=1;fragmts=1;match_lengths=298;repeat_ids=AluY;matching_classes=SINE/Alu;
      RM_hit_strands=+;RM_hit_IDs=72004;L1_5PINV=None;total_match_length=300;
      total_match_span=0.949367;ULTRA_TR=37;ULTRA_TR_span=0.117089;
      total_repeat_span=0.949367;polyA=TRUE
      GT  1|0  0/0  1|0
```

| Column | What GraffiTE puts there | Source |
|---|---|---|
| `ID` column | The caller's ID, with `_<n>` appended after the truvari merge (`n` is the record's rank in `SVs.vcf`): `chr1-18082-INS-315_10` above is PAV's `chr1-18082-INS-315`, tenth in the merge. svim-asm IDs are prefixed with the assembly name (`HG002_mat.svim_asm.INS.12`). With `--vcf` the IDs pass through unchanged. An ID longer than 50 characters stops Stage B. | <span class="src">`module/main.nf:154,232`, `bin/shorten_ids.py:21`, `bin/repmask_vcf.sh:13-16`</span> |
| `REF` and `ALT` columns | Sequence-resolved. For an insertion `REF` is the anchor base and `ALT` is that base plus the inserted sequence; for a deletion `REF` is the anchor base plus the deleted reference interval, re-read from the reference FASTA, and `ALT` is the anchor base. Symbolic `<INS>` and `<DEL>` records from Sniffles2 are dropped in Stage A. | <span class="src">`bin/fix_vcf.py:33-44`, `module/main.nf:98`</span> |
| `FILTER` column | Whatever the SV caller wrote. GraffiTE defines no FILTER of its own; see [Fields inherited from upstream callers](#fields-inherited-from-upstream-callers). | <span class="src">`module/main.nf:545`</span> |
| Sample columns | In `pangenome.vcf`, one column per Stage A sample, holding the genotype its caller wrote (PAV writes phased diploid calls, as above). A genotype missing after the merge is set to `0`. In `GraffiTE.merged.genotypes.vcf.gz`, one column per read set from `--genotype_with`. | <span class="src">`module/main.nf:230`, `module/main.nf:831`</span> |

Every VCF GraffiTE publishes carries a `##GraffiTE_version=` line right after `##fileformat`, with
the content of `version.txt` (`1.1.0`).
<span class="src">`module/main.nf:362,586,840`</span>

!!! warning "The ALT allele is not the element"
    For a `DEL` record the transposable element is in the reference and an `ALT` genotype means
    the sample lacks it. The presence-absence TSVs in `3_TSD_search/` fold this polarity for you;
    see [Output files](outputs.md).

## Repeat annotation fields

Added by `repeatmask_VCF` from the RepeatMasker and ULTRA output on each variant's inserted or
deleted sequence. A **hit** is one RepeatMasker element after its fragments have been grouped by
RepeatMasker's own link ID (column 15 of the `.out` file); a **fragment** is one line of that
file. Simple repeats and low-complexity fragments are removed before grouping, so they count
neither as hits nor toward the TE span.
<span class="src">`bin/annotate_vcf.R:73-90`</span>

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `n_hits` | 1 | Integer | Hits on the variant sequence. `0` when RepeatMasker found nothing. | <span class="src">`bin/repmask_vcf.sh:125`, `bin/annotate_vcf.R:148`</span> |
| `fragmts` | . | Integer | Fragments grouped into each hit, in hit order. | <span class="src">`bin/repmask_vcf.sh:129`, `bin/annotate_vcf.R:142`</span> |
| `match_lengths` | . | Integer | Bases of the variant covered by each hit, first to last fragment. | <span class="src">`bin/repmask_vcf.sh:126`, `bin/annotate_vcf.R:136`</span> |
| `repeat_ids` | . | String | Name of each hit, taken from its highest-scoring fragment. A hit whose fragments carry different names gets an `(x)` suffix. An SVA hit lying entirely inside the VNTR gets a `(VNTR_only)` suffix; see [SVA VNTR polymorphisms](../background/sva-vntr.md). | <span class="src">`bin/repmask_vcf.sh:127`, `bin/annotate_vcf.R:83-87,127-129`</span> |
| `matching_classes` | . | String | RepeatMasker class of each hit, as `class/family` from the library (`SINE/Alu`, `LINE/L1`, `LTR/ERVK`, ...). A `(VNTR_only)` SVA hit is reported as `Simple_repeat` here. | <span class="src">`bin/repmask_vcf.sh:128`, `bin/annotate_vcf.R:85,130-132`</span> |
| `RM_hit_strands` | . | String | Strand of each hit: `+` or `C` as RepeatMasker writes them. A hit whose fragments lie on both strands reports the concatenation (`C+`, `+C`), except an L1 with the `C+` twin-priming signature, whose strand is inferred; see [L1 5' inversions](../background/l1-5prime-inversion.md). | <span class="src">`bin/repmask_vcf.sh:130`, `bin/annotate_vcf.R:80,96-103`</span> |
| `RM_hit_IDs` | . | String | RepeatMasker link ID of each hit, so a hit can be found again in `repeatmasker_dir/indels.fa.out`. | <span class="src">`bin/repmask_vcf.sh:131`, `bin/annotate_vcf.R:146`</span> |
| `L1_5PINV` | . | String | Link IDs of the hits flagged as an L1 with a 5' inversion, or `None`. | <span class="src">`bin/repmask_vcf.sh:134`, `bin/annotate_vcf.R:95,147`</span> |
| `total_match_length` | 1 | Integer | Bases of the variant covered by TE hits, overlaps counted once. | <span class="src">`bin/repmask_vcf.sh:132`, `bin/repmask_vcf.sh:61-71`</span> |
| `total_match_span` | 1 | Float | `total_match_length` divided by the variant length. Written for continuity with v1.0, where it was the filter metric; nothing in v1.1 filters on it. | <span class="src">`bin/repmask_vcf.sh:133,71`</span> |
| `ULTRA_TR` | 1 | Integer | Bases of the variant that ULTRA annotates as tandem repeat, overlaps counted once. `0` when ULTRA found nothing. | <span class="src">`bin/repmask_vcf.sh:135,49-51`</span> |
| `ULTRA_TR_span` | 1 | Float | `ULTRA_TR` divided by the variant length, capped at 1. | <span class="src">`bin/repmask_vcf.sh:136,82-84`</span> |
| `total_repeat_span` | 1 | Float | Fraction of the variant covered by the union of TE hits and ULTRA intervals, capped at 1. This is the Stage B filter metric: records at or below `--repeat_span_cutoff` (default `0.80`) are discarded. | <span class="src">`bin/repmask_vcf.sh:137,86-94`, `module/main.nf:533,615`</span> |

A record without any hit has `n_hits=0`, `repeat_ids=None`, `matching_classes=None`,
`RM_hit_strands=None`, `RM_hit_IDs=None` and `L1_5PINV=None`; such records fail the span filter
unless ULTRA covers them.
<span class="src">`bin/annotate_vcf.R:164-171`</span>

## TSD and polyA fields

Added by `tsd_report` and `concat_repeatmask` to `pangenome.vcf`, and carried by every record
downstream, including the genotyped VCFs. A consolidated HERV-K record is built from scratch and
does not carry them.

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `TSD` | 1 | String | The target site duplication as its 5' copy and its 3' copy, comma-separated and upper-cased: `GATTACAG,GATTACAG`. The two copies are exact matches. Absent when the search found no duplication of 4 to 20 bp whose ends sit, on average, within 5 bp of the breakpoints. How the search works is in [Target site duplications](../background/tsd.md). | <span class="src">`bin/tsd_annotate_vcf.sh:18,27`</span> |
| `polyA` | 1 | String | `TRUE` when a tail of at least 8 bp with at least 80% A ends within 5 bp of the 3' end of a plus-strand hit, or the same with T at the 5' end of a minus-strand hit, after the matching `TSD` copy has been trimmed off. `FALSE` otherwise, including when the single hit has a mixed strand. `NA` when `n_hits` is greater than 1. | <span class="src">`bin/add_polyA.py:21-31,102-112`</span> |

## HERV-K classifier fields

Written only under `--human`, by `hervk_annotate`. `pangenome.vcf` is left untouched because it
induces the graph, so these fields appear in `pangenome.human.vcf`, in `hervk_candidates.vcf`
and in the two consolidated VCFs. What the states mean is in
[HERV-K (HML-2) biology](../background/hervk-hml2.md).

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `HERVK_CLASS` | 1 | String | Polymorphism class: `null_solo`, `solo_prov`, `truncated_prov`, `null_prov`, `copy_number` or `other`. | <span class="src">`bin/hervk_classify.py:425`</span> |
| `HERVK_ALLELE_REF` | 1 | String | State of the REF allele: `null`, `solo`, `provirus`, `prov_xN` (a tandem array of N proviral units sharing an LTR at each junction), `partial`, or `.` when unresolved. | <span class="src">`bin/hervk_classify.py:428`</span> |
| `HERVK_ALLELE` | . | String | State of each ALT allele, in ALT order, same vocabulary. | <span class="src">`bin/hervk_classify.py:432`</span> |
| `HERVK_EVIDENCE` | 1 | String | What resolved the states: `ARCH_2LTR`, `ARCH_PERM`, `ARCH_INT_PERM`, `CNV_PERIOD`, `ARCH_SOLO`, `REF_ANNOT`, `DENOVO_LTR`, `UNRESOLVED` or `NON_HML2`. | <span class="src">`bin/hervk_classify.py:434`</span> |
| `HERVK_ARCH` | 1 | String | 5' to 3' architecture of the variant allele with consensus intervals, for example `LTR:575-968/INT:1-7536/LTR:1-574`. | <span class="src">`bin/hervk_classify.py:437`</span> |
| `HERVK_K` | 1 | Integer | Where the aligner broke the reference solo LTR. An alignment property that varies between haplotypes and callers; nothing should key on its value. | <span class="src">`bin/hervk_classify.py:440`</span> |
| `HERVK_J` | 1 | Integer | The same for the internal region, the `ARCH_INT_PERM` counterpart of `HERVK_K`. | <span class="src">`bin/hervk_classify.py:443`</span> |
| `HERVK_N_UNITS_REF` | 1 | Integer | Proviral units in the masked reference element, a junction LTR counted once. | <span class="src">`bin/hervk_classify.py:446`</span> |
| `HERVK_UNIT_BP` | 1 | Integer | Period of the reference array in bp: one internal region plus one LTR. A copy-number change moves a whole number of these. | <span class="src">`bin/hervk_classify.py:448`</span> |
| `HERVK_REF_STATE` | 1 | String | State of the masked reference window: `null`, `solo`, `provirus`, `partial` or `unknown`. | <span class="src">`bin/hervk_classify.py:451`</span> |
| `HERVK_LAMBDA` | 1 | Float | bp of HML-2 LTR sequence on the variant allele. | <span class="src">`bin/hervk_classify.py:453`</span> |
| `HERVK_NU` | 1 | Float | bp of HML-2 internal (`HERVK-int`) sequence on the variant allele. | <span class="src">`bin/hervk_classify.py:455`</span> |
| `HERVK_COV` | 1 | Float | Fraction of the variant allele that is HML-2 sequence. | <span class="src">`bin/hervk_classify.py:457`</span> |
| `HERVK_PMAP` | 1 | Float | Confidence in the class under the size model. Reporting only; it does not decide the class. | <span class="src">`bin/hervk_classify.py:459`</span> |
| `HERVK_NOTE` | . | String | Diagnostics: `REF_ARCH_CONFLICT` (architecture and masked reference imply different REF states), `CNV_UNITS:a->b` (units on REF and ALT), `UNIT_COUNT_ASSUMED` (the count came from the architecture, since the reference did not resolve into counted units). | <span class="src">`bin/hervk_classify.py:462`</span> |

## HERV-K locus flags

Added by `hervk_reconcile.py flag` to `pangenome.human.vcf` after classification, and carried
onto the genotyped human VCF. They group records that describe the same element and say what
varies there. Records are never merged at this stage; the graph is induced from this file, so its
structure has to stay.

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `HERVK_LOCUS` | 1 | String | Locus identifier shared by the records that describe one element. | <span class="src">`bin/hervk_reconcile.py:96`</span> |
| `HERVK_LOCUS_N` | 1 | Integer | Records assigned to the locus. | <span class="src">`bin/hervk_reconcile.py:98`</span> |
| `HERVK_MEI` | 0 | Flag | A null allele segregates: the element is absent from some haplotypes, so a transposition produced the difference. This is the flag that separates an insertion polymorphism, comparable to an *Alu*, L1 or SVA insertion, from structural variation in an element every haplotype carries. | <span class="src">`bin/hervk_reconcile.py:100`</span> |
| `HERVK_SOLO_PROV` | 0 | Flag | A solo LTR and a provirus both segregate. Can be set beside the other two flags. | <span class="src">`bin/hervk_reconcile.py:106`</span> |
| `HERVK_CNV` | 0 | Flag | Some allele carries two or more proviral units. Can be set beside the other two flags. | <span class="src">`bin/hervk_reconcile.py:109`</span> |
| `HERVK_LOCUS_TYPE` | 1 | String | One label derived from the flags: `null_vs_present` when `HERVK_MEI` is set, else `copy_number`, else `solo_vs_provirus`, else `unresolved`. A locus that is two things at once keeps only the first in this field; the flags are the precise statement. | <span class="src">`bin/hervk_reconcile.py:113`</span> |
| `HERVK_MERGE_FLAG` | 0 | Flag | The locus holds more than one record. | <span class="src">`bin/hervk_reconcile.py:118`</span> |
| `HERVK_POLARITY_CONFLICT` | 0 | Flag | Records at the locus imply different REF states. The masked reference wins, and consolidation resolves it. | <span class="src">`bin/hervk_reconcile.py:122`</span> |

In the CaG cohort chr6:78,894,316 segregates a solo LTR, a provirus and a two-unit allele. It is
labelled `copy_number` and carries both `HERVK_SOLO_PROV` and `HERVK_CNV`.

## Consolidated HERV-K records

GraffiTE can emit several records for one HERV-K locus: two breakpoints for the same insertion,
or a deletion and an insertion describing opposite directions of the same event. The two
`.consolidated.` and `.human.` VCFs collapse each locus onto one multi-allelic record and add the
fields below. A `##hervk_consolidation=source:...` header line says whether the genotypes came
from discovery or from the graph.
<span class="src">`bin/hervk_reconcile.py:1164`</span>

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `SVTYPE` | A | String | Redefined per ALT, since one consolidated record can hold an insertion and a deletion. | <span class="src">`bin/hervk_reconcile.py:528`</span> |
| `SVLEN` | A | Integer | Redefined per ALT for the same reason. | <span class="src">`bin/hervk_reconcile.py:529`</span> |
| `HERVK_MEMBERS` | . | String | Record IDs consolidated into the locus. | <span class="src">`bin/hervk_reconcile.py:535`</span> |
| `HERVK_AC` | . | Integer | Allele count per ALT after dosage resolution across the members. | <span class="src">`bin/hervk_reconcile.py:537`</span> |
| `HERVK_AN` | 1 | Integer | Alleles called at the locus. Below `2N` it means members were structurally uncalled for some samples, since the graph can find every carrier and still fail to confirm the non-carriers. | <span class="src">`bin/hervk_reconcile.py:539`</span> |
| `HERVK_N_RESOLVED` | 1 | Integer | Samples fully resolved by dosage. | <span class="src">`bin/hervk_reconcile.py:542`</span> |
| `HERVK_N_PARTIAL` | 1 | Integer | Samples with some haplotypes unaccounted for; those haplotypes are reported missing. | <span class="src">`bin/hervk_reconcile.py:544`</span> |
| `HERVK_N_PLOIDY_EXCEEDED` | 1 | Integer | Samples whose summed ALT dosage exceeds their ploidy, which happens when a third allele was lost to `bcftools norm -m-`. Their genotypes are set missing. | <span class="src">`bin/hervk_reconcile.py:546`</span> |
| `HERVK_AC_DISC` | . | Integer | Allele count per ALT from the discovery (assembly) callset, an independent check on `HERVK_AC`. | <span class="src">`bin/hervk_reconcile.py:549`</span> |
| `HERVK_AN_DISC` | 1 | Integer | Alleles in the discovery callset. | <span class="src">`bin/hervk_reconcile.py:551`</span> |
| `HERVK_DISC_CONCORDANT` | 0 | Flag | Graph and discovery agree on every ALT count. | <span class="src">`bin/hervk_reconcile.py:553`</span> |
| `HERVK_GT_MASKED` | 0 | Flag | Graph genotypes were withheld at a copy-number locus (see the warning below). | <span class="src">`bin/hervk_reconcile.py:555`</span> |
| `HERVK_ALLELE_NOGT` | . | String | Alleles the graph cannot genotype, for which `HERVK_AC` reports 0. Their counts are in `HERVK_AC_DISC`. | <span class="src">`bin/hervk_reconcile.py:561`</span> |
| `HERVK_DISC_PLOIDY_MISMATCH` | 0 | Flag | Discovery and graph disagree on ploidy, so the discovery counts are withheld rather than reported on a denominator the two do not share. Expected on hemizygous chromosomes and with a haploid discovery caller such as svim-asm run per haplotype. | <span class="src">`bin/hervk_reconcile.py:567`</span> |
| `HERVK_MEMBERS_MASKED` | . | String | Members whose genotypes were withheld. | <span class="src">`bin/hervk_reconcile.py:572`</span> |
| `HERVK_ALLELE_SET` | . | String | Every allele state segregating at the locus, including states carried by members that are not in this file. | <span class="src">`bin/hervk_reconcile.py:574`</span> |
| `HERVK_MEMBERS_ABSENT` | . | String | Members the `--human` filter removed, whose alleles cannot be counted here. | <span class="src">`bin/hervk_reconcile.py:578`</span> |
| `HERVK_LOCUS_INCOMPLETE` | 0 | Flag | Some member is absent from this file, so the frequencies do not cover every allele. | <span class="src">`bin/hervk_reconcile.py:581`</span> |
| `HERVK_MEMBERS_UNRESOLVED` | . | String | Members with no usable allele state. | <span class="src">`bin/hervk_reconcile.py:584`</span> |
| `HERVK_POLARITY_FLIPPED` | . | String | Members whose own polarity disagreed with the locus REF state and were re-expressed against it. | <span class="src">`bin/hervk_reconcile.py:586`</span> |

`HERVK_LOCUS`, `HERVK_ALLELE_REF`, `HERVK_ALLELE` and the four locus flags are written again on
consolidated records with the same meaning.
<span class="src">`bin/hervk_reconcile.py:506-511,530-533`</span>

!!! warning "A copy-number allele cannot be genotyped from short reads"
    Its ALT path repeats sequence the REF path already carries, so a read from the pre-existing
    copy traverses it and non-carriers pick up ALT support. More depth does not remove that.
    GraffiTE withholds those alleles from the graph genotypes and names them in
    `HERVK_ALLELE_NOGT`. The other alleles at the same locus are genotyped normally, and the
    withheld one takes its frequency from the assemblies. At chr6:78,894,316 the graph counts
    the solo allele at AC=8, matching the assemblies, and the two-unit allele comes from
    `HERVK_AC_DISC`.

!!! warning "`HERVK_LOCUS_INCOMPLETE` means the frequencies are partial"
    The `--human` filter is narrower than the HERV-K candidate list, so it can remove some
    members of a locus and keep others. Three CaG loci are in that state. At chr7:4,699,714 it
    keeps one of three; all three describe the same 8.5 kb unit, and the two it removes carry a
    third RepeatMasker fragment the HERV-K clause does not admit. The surviving record carries
    the locus id and the full `HERVK_ALLELE_SET`, but its counts cover only the alleles still
    present. `hervk_candidates.vcf` has all the records with their assembly genotypes.

Three habits when reading these records. Quote `HERVK_AC` and `HERVK_AN` together rather than a
frequency alone, because `AN` below `2N` deflates the denominator without touching `AC`. Treat
`HERVK_N_PLOIDY_EXCEEDED` above zero as a sign of a lost third allele. Do not key on `HERVK_K`.

## FORMAT fields

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `GT` | 1 | String | Genotype. In `pangenome.vcf` it is the discovery genotype: `1` or `0` per haploid assembly, or the diploid call from Sniffles2. In the genotyped VCFs it is what PanGenie or `vg call` wrote for each read set. | <span class="src">`bin/repmask_vcf.sh:138`, `bin/merge_vcfs.py:521`</span> |

PanGenie and `vg call` add their own FORMAT fields (genotype quality, depth, likelihoods) with
their own header lines; GraffiTE passes them through unchanged. Which ones appear depends on the
version of those tools in the container, so they are not listed here.

## Fields inherited from upstream callers

`SVTYPE`, `SVLEN` and `END` come from the SV caller. With more than one input VCF every upstream
INFO field is stripped before the truvari merge, and `truvari collapse` adds `NumCollapsed`,
`NumConsolidated` and `CollapseId`; with a single input VCF the caller's INFO fields (PAV's `HAP`,
`QRY_REGION` and so on in the record above) are kept. `SVLEN` is then recomputed as
`strlen(ALT) - strlen(REF)` for every record, so it is negative for a deletion, and `abs(SVLEN)`
is what the trusted and human filters test.
<span class="src">`module/main.nf:186-196,231,482,498`</span>

GraffiTE defines no `##FILTER` line and never sets `FILTER` itself. `pangenome.vcf` retains the
value the caller wrote, and the caller's `##FILTER` definitions (PAV's `TRIM`, `COMPOUND`,
`QRY_FILTER` and `SVLEN`, for instance) travel with it. The trusted and human subsets require
`PASS` unless `--trusted_ignore_filter` or `--human_ignore_filter` is set.
<span class="src">`module/main.nf:483,527,545`</span>

One field is written to a file that is never published:

| Field | Number | Type | Meaning | Source |
|---|---|---|---|---|
| `ID` | A | String | Variant IDs per ALT allele, in the multi-allelic graph VCF that `pangenie_index` builds for PanGenie. Internal to the PanGenie path. | <span class="src">`bin/merge_vcfs.py:520`</span> |

## Retired fields

| Field (v1.0) | What replaced it |
|---|---|
| `mam_filter_1` (`5P_INV`) | `L1_5PINV`, written for every run. `--mammal` is accepted and ignored. |
| `mam_filter_2` (`VNTR_ONLY:...`) | The `(VNTR_only)` suffix on `repeat_ids` and `Simple_repeat` in `matching_classes`. |
| `total_match_span` as the filter metric | `total_repeat_span` with `--repeat_span_cutoff`. `total_match_span` is still written. |

<span class="src">`bin/repmask_vcf.sh:111-116`, `bin/annotate_vcf.R:127-132`</span>

## Filtering GraffiTE VCFs with bcftools

The multi-valued fields (`repeat_ids`, `matching_classes`, `RM_hit_strands`, everything with
`Number=.`) behave in two ways worth knowing before writing an expression. `~` is evaluated on
each element, and `^` anchors each element, so `repeat_ids~"^LTR5_Hs"` is true for a record
whose value is `SVA_A,LTR5_Hs`. And `!~` does not negate reliably on these fields:
`matching_classes!~"SINE"` is true for every record, `SINE/Alu` ones included. GraffiTE's own
`--human` expression is written with positive matches only for that reason. bcftools regular
expressions also have no alternation, which is why the whitelist parameters are comma-separated
lists rather than `A|B`.
<span class="src">`module/main.nf:486-503`</span>

Single-hit *Alu* insertions with a polyA tail and a TSD:

```bash
bcftools view -i 'n_hits=1 & matching_classes="SINE/Alu" & polyA="TRUE" & TSD!="."' pangenome.vcf
```

Read HERV-K from a consolidated VCF; in the unmerged ones a locus is several records and the same
allele can be counted twice. Insertion polymorphisms only, the set to use alongside *Alu*, L1 and
SVA:

```bash
bcftools view -i 'INFO/HERVK_MEI=1' GraffiTE.merged.genotypes.human.vcf.gz
```

Every other HERV-K locus, where the element is in every haplotype and its structure varies:

```bash
bcftools view -i 'INFO/HERVK_LOCUS!="." && INFO/HERVK_MEI=0' \
    GraffiTE.merged.genotypes.human.vcf.gz
```

Drop loci whose counts do not cover every allele, either because the graph could not count one or
because the `--human` filter removed a member:

```bash
bcftools view -e 'INFO/HERVK_ALLELE_NOGT!="." || INFO/HERVK_LOCUS_INCOMPLETE=1' \
    GraffiTE.merged.genotypes.human.vcf.gz
```

SVA VNTR-only records, which carry `Simple_repeat` rather than `Retroposon/SVA`:

```bash
bcftools view -e 'INFO/repeat_ids~"VNTR_only"' pangenome.human.vcf
```
