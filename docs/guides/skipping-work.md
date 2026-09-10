---
title: Resuming and skipping work
description: >-
  Nextflow -resume, and the eight entry points that let a run start from files an earlier run
  produced instead of recomputing them.
---

# Resuming and skipping work

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

Two mechanisms avoid recomputing. Nextflow's `-resume` reuses cached tasks from the same launch
directory when nothing upstream changed. GraffiTE's own entry-point flags go further: they take
a file an earlier run published, or that you produced elsewhere, and skip every process before
it. The flags form a ladder, from the earliest stage to the latest.

---

## Nextflow `-resume`

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -latest -resume --assemblies assemblies.csv ...
```

Nextflow keys each task on its inputs and script, so a task reruns only when one of those
changed. This works only from the directory holding the `work/` tree and `.nextflow/` of the
earlier run, and only as long as `work/` has not been deleted. It is the right tool for a run
that stopped part-way. The flags below are for reusing results across directories, across
machines, or with different downstream parameters.

---

## The ladder

| Flag | What you supply | Skips | Source |
|---|---|---|---|
| `--svs` | per-sample SV VCFs | the callers, not the merge | <span class="src">`main.nf:120-125`</span> |
| `--vcf` | one merged SV VCF | all of Stage A | <span class="src">`main.nf:148-149`</span> |
| `--RM_dir` | the `2_Repeat_Filtering/` directory of an earlier run | RepeatMasker and ULTRA | <span class="src">`main.nf:133-142`</span> |
| `--graffite_vcf` | a `pangenome.vcf` | all of Stages A and B | <span class="src">`main.nf:129, 183-186`</span> |
| `--graph` | a `GraffiTE_graph/index/` directory | `make_graph` | <span class="src">`main.nf:221-225`</span> |
| `--graph_alignments` | per-sample GAF and pack files | `graph_align_reads` | <span class="src">`main.nf:234-240`</span> |
| `--vcfs` | per-sample `vg call` VCFs | alignment and `vg_call` | <span class="src">`main.nf:229-231`</span> |
| `--hervk_reconcile_vcf` | a genotyped VCF | genotyping, with `--genotype false` | <span class="src">`main.nf:292-310`</span> |
| `--genotype false` | nothing | all of Stage C | <span class="src">`main.nf:188`</span> |

Each rung is described below with what must exist on disk and what still runs.

### `--svs`: your own per-sample calls

**On disk:** a samplesheet with `sample` and `path` columns, one VCF per row; see
[Samplesheets](../reference/samplesheets.md).

**Runs:** the merge (`truvari_merge`) and everything after it. Combines with `--assemblies`,
`--longreads`, `--bams` and `--pav`; the VCFs are mixed into the same merge. GraffiTE does not
filter what you supply by size or type, so restrict to `INS` and `DEL` yourself
<span class="src">`main.nf:120-125`</span>.

### `--vcf`: one merged VCF

**On disk:** one VCF or `.vcf.gz` of insertions and deletions, with sequence-resolved alleles.
Multi-allelic records should be split first (`bcftools norm -m-`).

**Runs:** `truvari_merge` in pass-through mode (decompress only, IDs preserved), then all of
Stage B and C <span class="src">`module/main.nf:171-178`</span>.

**Cannot be combined** with any discovery flag; the run stops with a message naming the flag it
saw <span class="src">`main.nf:54-57`</span>. Use `--svs` when you want your VCF merged with
GraffiTE's own calls.

### `--RM_dir`: RepeatMasker output of an earlier run

**On disk:** a directory holding one subdirectory per Stage B batch, each with
`genotypes_repmasked_filtered.vcf` and a `repeatmasker_dir/`. That is exactly the layout of
`out/2_Repeat_Filtering/`, where the subdirectories are numbered by task
<span class="src">`main.nf:134-136`, `module/main.nf:592-598`</span>.

```text
out/2_Repeat_Filtering/
├── 1/
│   ├── genotypes_repmasked_filtered.vcf
│   └── repeatmasker_dir/
├── 2/
│   ...
```

**Runs:** the TSD search, polyA annotation, the repeat-span filter, the trusted or human
subset, and, with `--human`, the HERV-K classifier, which reads the raw RepeatMasker tables from
`repeatmasker_dir/`. `--TE_library` is still required under `--human`, because the classifier
masks reference windows with it <span class="src">`main.nf:172-175`</span>.

This is the entry point for changing anything downstream of RepeatMasker: `--tsd_win`,
`--repeat_span_cutoff`, the `--trusted_*` and `--human*` parameters, or the HERV-K parameters.

### `--graffite_vcf`: a finished `pangenome.vcf`

**On disk:** a `pangenome.vcf` from `3_TSD_search/`, or any VCF carrying the GraffiTE INFO
fields.

**Runs:** Stage C only, which builds the graph from the VCF as-is
<span class="src">`main.nf:183-186`</span>.

**Cannot be combined** with `--human` unless `--hervk_reconcile false` is also given: the HERV-K
consolidation after genotyping needs tables that only the skipped stage writes
<span class="src">`main.nf:66-71`</span>. To re-run the HERV-K steps, enter with `--RM_dir`
instead.

### `--graph`: a built graph

**On disk:** the `GraffiTE_graph/index/` directory of an earlier run, or any directory holding
what `make_graph` writes: `index.gfa` and `index.pb`, plus `index.giraffe.gbz` and its
companions for the `giraffe` method <span class="src">`module/main.nf:723-738`</span>.

**Runs:** alignment and calling for every sample in `--genotype_with`. Requires
`--graph_method giraffe`, `graphaligner` or `precomputed`; the PanGenie path has its own index
and ignores `--graph` <span class="src">`main.nf:214-225`</span>.

### `--graph_alignments`: aligned reads

**On disk:** a samplesheet with `sample`, `gaf` and `pack` columns pointing at the
`<sample>.gaf.gz` and `<sample>.pack` files from `GraffiTE_alignments/`
<span class="src">`main.nf:234-236`</span>.

**Runs:** `vg_call` for each row, against the graph from `--graph` or a fresh `make_graph`.
The rows here decide which samples are called; `--genotype_with` is still read, and every path
in it must exist, but its samples are not aligned.

### `--vcfs`: called samples

**On disk:** a samplesheet with `sample` and `path` columns pointing at per-sample `vg call`
VCFs (`.vcf.gz` with index), as `vg_call` writes them <span class="src">`main.nf:229-231`</span>.

**Runs:** `merge_VCFs` only, and `hervk_reconcile` under `--human`. No sample is aligned or
called, though `--genotype_with` is still read and its paths must exist; `--graph` is still needed for the `precomputed` method's validation
<span class="src">`main.nf:59-64`</span>.

### `--hervk_reconcile_vcf` with `--genotype false`

**On disk:** a `GraffiTE.merged.genotypes.vcf.gz` from an earlier run, and the inputs for
`--RM_dir` and `--human`.

**Runs:** Stage B from the RepeatMasker output, `hervk_annotate`, then `hervk_reconcile`
against the supplied VCF, whose back end is read from its header
<span class="src">`main.nf:298-310`</span>. See [Human MEIs](human-mei.md) for the command
line.

### `--genotype false`

Stops after Stage B. `--genotype_with` is not read and need not exist
<span class="src">`main.nf:188`</span>.

---

## Combining rungs

Rungs at different stages combine: `--RM_dir` with `--graph`, or `--graffite_vcf` with
`--graph_alignments`. The rule is that each flag replaces the processes that would have produced
its input and leaves the rest running, so any combination that leaves a connected pipeline works.
The one refusal is `--vcf` beside a discovery flag.

Rungs at the same stage do not stack: `--graph_alignments` and `--vcfs` together means `--vcfs`
wins, because the alignment branch is inside the `else` of the `--vcfs` test
<span class="src">`main.nf:229-243`</span>.
