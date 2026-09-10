---
title: Human mobile element insertions
description: >-
  The --human polymorphic MEI filter, gate by gate, and the HERV-K (HML-2) allele-state
  classifier and locus consolidation that run behind it.
---

# Human mobile element insertions

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

`--human` does two things. In Stage B it writes `pangenome.human.vcf`, a subset of
`pangenome.vcf` restricted to the mobile element subfamilies still active in humans, **instead
of** `pangenome.trusted.vcf` <span class="src">`module/main.nf:551-562`</span>. Then it runs the
HERV-K (HML-2) classifier over the callset, and after genotyping it consolidates the HERV-K loci
in the genotyped calls.

The filter reads `pangenome.vcf` directly; none of the `--trusted_*` parameters affect it. Class-level filtering would not be enough: `SINE/Alu` alone admits AluJ and AluS elements
that have not been mobile for tens of millions of years, so every class carries a subfamily
whitelist.

---

## What `--human` changes

<figure>
--8<-- "assets/human-funnel.svg"
<figcaption>
Every record of <code>pangenome.vcf</code> passes five gates in one <code>bcftools view -i</code>
expression. The carve-out on the right admits HML-2 proviral records that RepeatMasker split
into more than one hit.
</figcaption>
</figure>

The whole filter is one expression, assembled in Groovy and written verbatim into
`3_TSD_search/human_filter_summary.txt` together with record counts and the kept and dropped
`(matching_classes, repeat_ids)` combinations <span class="src">`module/main.nf:564-579`</span>.
Read that file first when a record you expected is missing.

---

## The subfamily whitelists

A record of a class passes only when one of its `repeat_ids` matches that class's whitelist regex. The
whitelists are comma-separated lists of bcftools regexes; each list is expanded into an `OR` of
`repeat_ids~"<regex>"` clauses <span class="src">`module/main.nf:491-497`</span>.

| `matching_classes` | Parameter | Default | Keeps | Source |
|---|---|---|---|---|
| `SINE/Alu` | `--human_alu_ids` | `^AluY` | AluY and its subfamilies (AluYa5, AluYb8, ...); drops AluS and AluJ | <span class="src">`nextflow.config:65`</span> |
| `LINE/L1` | `--human_l1_ids` | `^L1HS` | L1HS only; add `^L1PA2` to relax | <span class="src">`nextflow.config:66`</span> |
| `Retroposon/SVA` | `--human_sva_ids` | `^SVA_[DEF]` | SVA_D, SVA_E, SVA_F | <span class="src">`nextflow.config:67`</span> |
| `Simple_repeat` | `--human_sva_ids` | `^SVA_[DEF]` | the `SVA_*(VNTR_only)` records that Stage B reclassified as `Simple_repeat`; see [SVA VNTR polymorphisms](../background/sva-vntr.md) | <span class="src">`module/main.nf:496`</span> |
| `LTR/ERVK` | `--human_hervk_ids` | `^HERVK-int,^HERVK$,^LTR5_Hs,^LTR5A,^LTR5B` | the HML-2 internal region under either library name, and the LTR5 family; `^HERVK$` is anchored at both ends so that HERVK9-int, HERVK11-int and HERVK14-int, which are other lineages, stay out | <span class="src">`nextflow.config:68-76`</span> |

An empty string keeps the whole class <span class="src">`module/main.nf:492`</span>.

!!! warning "Two bcftools behaviours the filter depends on"
    `repeat_ids` and `matching_classes` are `Number=.` fields. bcftools evaluates `~` on them
    **element by element**, with `^` anchoring each element, which is why `repeat_ids~"^LTR5_Hs"`
    matches the record `SVA_A,LTR5_Hs`. bcftools regexes have **no alternation**, so each
    whitelist is a list rather than one `a|b` pattern. Negation with `!~` is not reliable on these
    fields either, which is why every clause below is phrased positively
    <span class="src">`module/main.nf:485-503`</span>. Details in
    [VCF fields](../reference/vcf-fields.md).

---

## Length and tandem-repeat gates

<span class="src">`module/main.nf:498`</span>

```text
abs(SVLEN) >= 250  and  (ULTRA_TR_span < 0.6  or  matching_classes = "Simple_repeat")
```

| Parameter | Default | Effect |
|---|---|---|
| `--human_min_svlen` | `250` bp | drops records shorter than a truncated Alu |
| `--human_max_ultra_span` | `0.6` (fraction of the variant covered by ULTRA tandem repeats) | drops VNTR-dominated records, except the `Simple_repeat` ones that are VNTR by construction |

Source: <span class="src">`nextflow.config:77-78`</span>.

---

## One element, or the HERVK + SVA carve-out

<span class="src">`module/main.nf:504-525`</span>. A record passes the structure gate in one of
two ways.

**One hit.** `n_hits == 1`, and the record is `LTR/ERVK`, `Simple_repeat`, or carries
`polyA = TRUE`. Alu, L1 and SVA insertions arrive by target-primed reverse transcription and
carry a polyA tail; HML-2 elements and VNTR expansions do not, so the filter requires the tail
only where the mechanism predicts it.

**The carve-out.** RepeatMasker often reports a small `SVA_*` hit at the end of an HML-2
provirus. That is homology, not a second element: SVA's SINE-R domain derives from the HERV-K
LTR (see [HERV-K biology](../background/hervk-hml2.md)). Such records have two hits and would
fail the one-hit rule, so the carve-out admits them when all of these hold:

```text
n_hits <= 3  and  LTR/ERVK  and  Retroposon/SVA  and  repeat_ids ~ "^HERVK-int"  and  abs(SVLEN) <= 10500
```

| Parameter | Default | Source |
|---|---|---|
| `--hervk_sva_pair` | `true`; `false` removes the carve-out | <span class="src">`nextflow.config:80`</span> |
| `--hervk_pair_max_svlen` | `10500` bp, a 9,472 bp provirus plus tolerance | <span class="src">`nextflow.config:81`</span> |
| `--hervk_pair_max_hits` | `3`; `2` restores the pre-1.1 rule | <span class="src">`nextflow.config:82`</span> |

The hit count is a cap rather than an equality because RepeatMasker also splits the internal
region of a degraded or rearranged provirus. In the CaG cohort, two of the three records at
chr7:4,699,714 came back with three hits (`LTR5_Hs,SVA_A,HERVK-int` and
`HERVK-int,SVA_A,HERVK-int`), where the third hit is that split. At `n_hits == 2` the filter
kept one record of the locus and dropped the two carrying its common allele, so the locus
reported a two-unit against three-unit difference where the assemblies show three states. The
other three clauses keep the cap honest: relaxing `n_hits` on its own in the one-hit rule admits
records that are an LTR5 fragment beside something else, alpha satellite and HERVK9 among them
<span class="src">`module/main.nf:509-523`</span>.

!!! note "The carve-out spells the internal region `HERVK-int`"
    The clause matches `repeat_ids~"^HERVK-int"` literally, not the `--human_hervk_ids` list. A
    library that names the internal region `HERVK`, as Dfam does, satisfies the whitelist but not
    the carve-out, so its split proviral records fall back on the one-hit rule.

Finally, unless `--human_ignore_filter` is set, the record must carry `FILTER = PASS` from the
SV caller <span class="src">`module/main.nf:527`</span>.

---

## The HERV-K classifier

Alu, L1 and SVA records are done once they pass the filter. HML-2 records are not: a solo LTR
and a full provirus at the same site are two different alleles of one locus, and a structural
variant call only ever describes the difference between two of a locus's states. `hervk_annotate`
works out which two, for every `LTR/ERVK` record in `pangenome.vcf`, and writes the result into
`pangenome.human.vcf` <span class="src">`main.nf:171-182`, `module/main.nf:252-366`</span>. The
biology is on the [HERV-K background page](../background/hervk-hml2.md); this is the procedure.

1. **Candidates.** Every record with `LTR/ERVK` among its classes and `|SVLEN|` at most
   `--hervk_max_svlen` (`25000` bp). This runs over the whole of `pangenome.vcf`, not only the
   human subset, so that a locus split by the filter is still seen whole
   <span class="src">`module/main.nf:300-301`</span>.
2. **Architecture.** `hervk_arch.py` re-reads the raw RepeatMasker tables (the `repeat_ids`
   field has collapsed each link group to one name, which erases the LTR-INT-LTR order), tiles
   the hits along the SV allele, reassigns SINE-R hits to the LTR class, and reports a signature:
   `ARCH_2LTR` (two complete terminal LTRs), `ARCH_PERM` (one LTR split across the two ends at
   permutation point `k`), `ARCH_INT_PERM` (the same split inside the internal region, point
   `j`), `ARCH_SOLO` (one LTR, no internal region) or `ARCH_NONE`. Written to `hervk_arch.tsv`
   <span class="src">`module/main.nf:303`, `bin/hervk_arch.py:254-356`</span>.
3. **Reference state.** `hervk_ref_state.py` cuts a window of the reference around each candidate
   (`--hervk_ref_flank`, `1500` bp each side), masks it with RepeatMasker against the TE library
   on `task.cpus` threads, and reads off what the reference holds: `null`, `solo`, `provirus`,
   `partial`, or `unknown`. An element that runs into the window edge is re-masked with a wider
   flank, up to three doublings of 12,000 bp. For a proviral reference it counts the units and
   measures their period. With `--hervk_ref_annotation` (a RepeatMasker `.out` or a BED of the
   reference) the classifier skips the masking step. Written to `hervk_refstate.tsv`
   <span class="src">`module/main.nf:305-314`, `bin/hervk_ref_state.py:46-72`</span>.
4. **Classification.** `hervk_classify.py` combines the two tables and `SVLEN` into
   `HERVK_ALLELE_REF`, `HERVK_ALLELE`, `HERVK_EVIDENCE` and `HERVK_CLASS`. Copy-number
   arithmetic runs first: when `|SVLEN|` is a whole number of the reference element's period and
   one side has two or more units, the states are `prov_xN`. Otherwise the architecture settles
   the pair, and where it is degenerate the reference state decides. Size alone never decides a
   class; `HERVK_PMAP` is a confidence under a size model, reported for information. It runs twice:
   over the full candidate set, writing `hervk_candidates.vcf`, `hervk_calls.tsv` and
   `hervk_polymorphism_summary.md`; and over the human subset, annotating its records and its
   presence-absence TSV <span class="src">`module/main.nf:323-334`, `bin/hervk_classify.py:290-391`</span>.
5. **Loci.** `hervk_reconcile.py flag` groups records whose footprints lie within
   `--hervk_locus_window` (`1200` bp, one LTR plus tolerance) into loci and writes
   `HERVK_LOCUS`, `HERVK_LOCUS_N`, the flags `HERVK_MEI`, `HERVK_SOLO_PROV`, `HERVK_CNV`,
   `HERVK_MERGE_FLAG`, `HERVK_POLARITY_CONFLICT`, and the summary `HERVK_LOCUS_TYPE`. The flag step marks records and never
   merges them, because `pangenome.human.vcf` must keep its record structure. The locus
   table is `hervk_loci.tsv` <span class="src">`module/main.nf:336-340`, `bin/hervk_reconcile.py:95-170`</span>.
6. **Consolidated discovery VCF.** The same loci collapsed to one multi-allelic record each, with
   the allele set and the counts from the assemblies, as a separate file
   `pangenome.human.consolidated.vcf` and a report <span class="src">`module/main.nf:353-359`</span>.

`pangenome.vcf` itself is read and never rewritten here: it induces the graph and must stay
byte-identical <span class="src">`module/main.nf:247-248`</span>. `pangenome.human.vcf` and its
TSV, first written by `concat_repeatmask`, are overwritten by this process's annotated versions
<span class="src">`module/main.nf:253`</span>.

The field-by-field definitions are in [VCF fields](../reference/vcf-fields.md).

---

## After genotyping: `hervk_reconcile`

Genotyping does not preserve the HERV-K annotation, because `merge_VCFs` copies INFO from the
un-annotated `pangenome.vcf`. `hervk_reconcile` restores it and consolidates the loci in the
genotyped calls <span class="src">`main.nf:282-289`, `module/main.nf:374-445`</span>:

1. `GraffiTE.merged.genotypes.vcf.gz` is subset to the IDs in `pangenome.human.vcf`.
2. Every `HERVK_*` INFO field is copied across from the discovery VCF with `bcftools annotate`.
3. `hervk_reconcile.py consolidate` writes one record per locus, with `HERVK_AC`/`HERVK_AN`
   from the graph genotypes beside `HERVK_AC_DISC`/`HERVK_AN_DISC` from the assemblies, and
   withholds the graph genotypes at copy-number loci (`HERVK_GT_MASKED`). At those loci the ALT
   path repeats sequence the reference already carries, so reads from the pre-existing copy
   traverse it and non-carriers acquire ALT support; the discovery genotypes come from
   haplotype-resolved alignments and do not have this problem
   <span class="src">`nextflow.config:98-106`</span>.

Outputs: `4_Genotyping/GraffiTE.merged.genotypes.human.vcf.gz` (indexed),
`hervk_unconsolidated_records.vcf` (the member records, archived) and
`hervk_reconciliation_report.md`. The full `GraffiTE.merged.genotypes.vcf.gz` is never rewritten.

!!! warning "Giraffe is the only validated back end"
    The reconciler refuses any other genotyper rather than produce an unchecked answer
    <span class="src">`bin/hervk_reconcile.py:431, 769-773`</span>. With `--graph_method
    pangenie` or `graphaligner`, `hervk_reconcile` exits with an error at the end of an otherwise
    complete run; pass `--hervk_reconcile false` to skip it. When the VCF comes in through
    `--hervk_reconcile_vcf`, the back end is detected from the header, and any `vg call` header is
    reported as `giraffe`.

### Re-running the HERV-K stages without genotyping

Graph genotyping is the expensive stage, and the HERV-K work does not depend on how the calls
were made. To repeat stages 1 to 6 and the consolidation against a genotyped VCF from an earlier
run <span class="src">`main.nf:292-310`</span>:

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -latest \
  --RM_dir earlier/out/2_Repeat_Filtering \
  --reference ref.fa --TE_library TEs.fa --human \
  --genotype false \
  --hervk_reconcile_vcf earlier/out/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz
```

`--RM_dir` is the entry point to use here, not `--graffite_vcf`: the classifier needs the raw
RepeatMasker tables, and `--graffite_vcf` skips the stage that provides them. The workflow
refuses `--graffite_vcf --human` unless `--hervk_reconcile false` is also given
<span class="src">`main.nf:66-71`</span>.

---

## Tuning

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--hervk_config` | `null` | JSON file of overrides for `hervk_classify.py`; keys present replace the script's `DEFAULTS`, nested dictionaries are merged key by key | <span class="src">`nextflow.config:61`, `bin/hervk_classify.py:58-97`</span> |
| `--hervk_strict` | `false` | drop candidates classed `other` or below `pmap_min` (`0.90`); off because dropping records hid a failure in the previous classifier | <span class="src">`nextflow.config:115`</span> |
| `--hervk_max_svlen` | `25000` bp | candidacy cap; a complete provirus is 9,472 bp and a 25 Mb artefact once dominated the masking cost | <span class="src">`nextflow.config:84`</span> |
| `--hervk_ref_flank` | `1500` bp | reference masked either side of each candidate footprint | <span class="src">`nextflow.config:108`</span> |
| `--hervk_locus_window` | `1200` bp | maximum footprint gap within one locus | <span class="src">`nextflow.config:109`</span> |
| `--hervk_ref_annotation` | `null` | precomputed reference repeat track; skips the in-pipeline masking | <span class="src">`nextflow.config:113`</span> |
| `--hervk_mask_graph_gt_at_cnv` | `true` | withhold graph genotypes at copy-number loci; `--hervk_mask_tandem` is the deprecated name | <span class="src">`nextflow.config:98-107`</span> |
| `--hervk_reconcile` | `true` | run the consolidation after genotyping | <span class="src">`nextflow.config:94`</span> |
| `--hervk_annotate_threads`, `_memory`, `_time` | `1`, `10G`, `12h` | resources for the masking step | <span class="src">`nextflow.config:127-129`</span> |

The keys `--hervk_config` can override, with their defaults, are the `DEFAULTS` dictionary at
the top of `bin/hervk_classify.py`: `sigmas`, `priors`, `t_min`, `t_max`, `int_full_frac`,
`pmap_min`, `min_hml2_bp`, `hml2_frac_min`, `cnv_period_tol`, `cnv_period_frac`, `max_svlen`,
`min_solo_bp` and `min_prov_int_bp`. Unknown keys are accepted and ignored.

!!! warning "`utils/HERVK.config.json` predates this classifier"
    The template in `utils/` carries keys from the previous, size-based classifier (`s_C`,
    priors named `C`/`T`/`B`/`A`/`X`, `sva_mimic_*`). None of them is read by the current
    script. Write your JSON against the `DEFAULTS` names above.

---

## Outputs specific to `--human`

All in `3_TSD_search/` unless stated; see [Output files](../reference/outputs.md) for the
directory tree and [VCF fields](../reference/vcf-fields.md) for the fields.

| File | Written by | Contents |
|---|---|---|
| `pangenome.human.vcf` | `concat_repeatmask`, then overwritten by `hervk_annotate` | the filtered subset, with HERV-K allele states and locus flags; this is what induces the graph together with `pangenome.vcf` |
| `pangenome.presence-absence_human.tsv` | same | its presence-absence table |
| `human_filter_summary.txt` | `concat_repeatmask` | the filter expression, counts, kept and dropped combinations |
| `hervk_candidates.vcf` | `hervk_annotate` | every HERV-K candidate in `pangenome.vcf`, annotated, whether or not the filter kept it |
| `hervk_calls.tsv`, `hervk_loci.tsv` | `hervk_annotate` | the per-record calls and the locus table |
| `hervk_arch.tsv`, `hervk_refstate.tsv` | `hervk_annotate` | the architecture and reference-state tables the calls rest on |
| `hervk_polymorphism_summary.md` | `hervk_annotate` | counts per class and evidence |
| `pangenome.human.consolidated.vcf`, `hervk_discovery_consolidation_report.md` | `hervk_annotate` | one record per locus, discovery genotypes |
| `4_Genotyping/GraffiTE.merged.genotypes.human.vcf.gz` | `hervk_reconcile` | the genotyped human subset, loci consolidated |
| `4_Genotyping/hervk_unconsolidated_records.vcf`, `hervk_reconciliation_report.md` | `hervk_reconcile` | the archived member records and the report |

`pangenome.trusted.vcf` and its TSV are **not** written under `--human`
<span class="src">`module/main.nf:551-555`</span>.
