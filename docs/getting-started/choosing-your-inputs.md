---
title: Choosing your inputs
description: >-
  Which GraffiTE entry flag matches the data you have, which stages each one
  runs, and which flags can be combined.
---

# Choosing your inputs

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `e687712`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](v1.0-vs-v1.1.md).

GraffiTE is three stages in series, and every input flag is a way of entering that series at a
given point. Choose the flag by asking what you already have.

---

## Stages and entry flags

| Stage | What it finds | Enters with | Produces |
|---|---|---|---|
| **1, SV search** | Where the sequence-resolved insertions and deletions sit relative to the reference | `--assemblies`, `--pav`, `--longreads`, `--bams`, `--svs` | `1_SV_search/SVs.vcf` |
| **2, TE annotation** | Which of those SVs are transposable elements, and what they are | `--vcf`, `--RM_dir` | `3_TSD_search/pangenome.vcf` and its subsets |
| **3.1, genotyping** | Which samples carry each polymorphism | `--graffite_vcf`, `--graph`, `--graph_alignments`, `--vcfs` | `4_Genotyping/GraffiTE.merged.genotypes.vcf.gz` |

The optional Stage 3.2 lifts methylation onto the graph (`--epigenomes`), only on the
`giraffe`, `graphaligner` and `precomputed` methods <span class="src">`main.nf:206,234`</span>. See
[Epigenetic mapping](../guides/methylation.md).

---

## What data do you have?

```mermaid
%%{init: {"flowchart": {"nodeSpacing": 14, "rankSpacing": 26, "padding": 8}}}%%
flowchart LR
    Q0("What do you<br/>start from?")
    Q0 -->|"Genome<br/>assemblies"| Q1("Call them<br/>with PAV?")
    Q1 -->|No| SA["<code>--assemblies</code>"]
    Q1 -->|Yes| PV["<code>--pav</code>"]
    Q0 -->|"Long reads"| Q2("Already aligned?")
    Q2 -->|No| LR["<code>--longreads</code>"]
    Q2 -->|Yes| BM["<code>--bams</code>"]
    Q0 -->|"SV calls with<br/>REF and ALT<br/>sequences"| Q3("One VCF per sample,<br/>or one merged VCF?")
    Q3 -->|"Per sample"| SV["<code>--svs</code>"]
    Q3 -->|"Merged"| VC["<code>--vcf</code>"]
    Q0 -->|"An earlier<br/>GraffiTE run"| Q4("What did<br/>it finish?")
    Q4 -->|"RepeatMasker"| RM["<code>--RM_dir</code>"]
    Q4 -->|"pangenome.vcf"| GV["<code>--graffite_vcf</code>"]
    Q4 -->|"Graph, alignments<br/>or vg call VCFs"| PC["<code>--graph</code><br/><code>--graph_alignments</code><br/><code>--vcfs</code>"]

    classDef q fill:#8888881a,stroke:#888888;
    classDef s1 fill:#f5851f2e,stroke:#f5851f;
    classDef s2 fill:#f2523f2e,stroke:#f2523f;
    classDef s3 fill:#3d63c92e,stroke:#3d63c9;
    class Q0,Q1,Q2,Q3,Q4 q;
    class SA,PV,LR,BM,SV s1;
    class VC,RM s2;
    class GV,PC s3;
```

<p class="gt-key"><span class="gt-key-s1">enters at Stage 1</span> <span class="gt-key-s2">at Stage 2</span> <span class="gt-key-s3">at Stage 3</span></p>

Every discovery input needs `--reference` and `--TE_library`. Genotyping needs `--genotype_with`
unless you pass `--genotype false`.

---

## Command builder

Answer the questions to get a starting command. Its file names are those of the examples on
[Samplesheet formats](../reference/samplesheets.md); replace them with your own paths.

<div class="gt-picker" data-revision="v1.1dev" markdown>
<p>The command builder needs JavaScript; without it, use the diagram above and the tables
below.</p>

[Samplesheet formats](../reference/samplesheets.md){ data-key="sheets" hidden }
[Resuming and skipping work](../guides/skipping-work.md){ data-key="skipping" hidden }
[Stage 1: SV search](../guides/discovery.md){ data-key="discovery" hidden }
[Human mobile element insertions](../guides/human-mei.md){ data-key="human" hidden }
[Stage 3.2: epigenetic mapping](../guides/methylation.md){ data-key="methylation" hidden }
[Installation](installation.md){ data-key="install" hidden }

</div>

---

## What each flag runs and skips

| Flag | Takes | Runs | Skips |
|---|---|---|---|
| `--assemblies` | CSV `sample,path` of haploid assemblies | minimap2 or winnowmap under `--aligner`, svim-asm, then Stages 2 and 3 | nothing |
| `--pav` | CSV read by position, first line skipped as a header, `sample,hap1[,hap2...]` | PAV in its own container, then Stages 2 and 3 | nothing |
| `--longreads` | CSV `sample,path,type` | minimap2 or winnowmap, Sniffles2, then Stages 2 and 3 | nothing |
| `--bams` | CSV `sample,path` of long-read alignments | Sniffles2, then Stages 2 and 3 | read alignment |
| `--svs` | CSV `sample,path` of per-sample VCFs | the truvari merge, then Stages 2 and 3 | all callers |
| `--vcf` | one sequence-resolved VCF | Stages 2 and 3 | the callers and the truvari collapse |
| `--RM_dir` | a `2_Repeat_Filtering` directory | TSD search onwards, then Stage 3 | Stage 1 and RepeatMasker |
| `--graffite_vcf` | a `pangenome.vcf` | Stage 3 only | Stages 1 and 2 |
| `--graph` | a `GraffiTE_graph/index` directory | read alignment and calling | graph construction |
| `--graph_alignments` | CSV `sample,gaf,pack` | `vg call` and the merge | graph alignment |
| `--vcfs` | CSV `sample,path` of per-sample `vg call` VCFs | the merge | alignment and calling |
| `--hervk_reconcile_vcf` | a genotyped VCF from an earlier run, with `--human --genotype false` | Stages 1 and 2, then HERV-K consolidation against that VCF | Stage 3 |
| `--genotype false` | | Stages 1 and 2 | Stage 3 |

Lines in `main.nf`: the discovery block <span class="src">`main.nf:76-126`</span>, the
annotation entry points <span class="src">`main.nf:129-156`</span>, `--graffite_vcf`
<span class="src">`main.nf:183-186`</span>, the genotyping block
<span class="src">`main.nf:188-294`</span>. Samplesheet columns are on
[Samplesheet formats](../reference/samplesheets.md).

---

## Combining flags

**Discovery flags add up.** You can pass `--assemblies`, `--pav`, `--longreads`, `--bams` and
`--svs` together in any combination; every caller's output goes into one truvari merge
<span class="src">`main.nf:125`</span>. The paper's `GT-svsn` mode is one such combination.

**`--vcf` stands alone.** It replaces the merge, so GraffiTE refuses any discovery flag beside it
before anything runs:

```
--vcf cannot be combined with --assemblies. Pass --vcf alone, or drop it and use --svs to add your own per-sample VCFs to the discovery merge.
```

<span class="src">`main.nf:46-49`</span>.

**`--graffite_vcf` skips the HERV-K annotation.** With `--human`, the consolidation step needs
files that only the annotation stage writes, so with `--graffite_vcf --human` pass
`--hervk_reconcile false`, or start from `--RM_dir` instead <span class="src">`main.nf:61-63`</span>.

**`--graph_method precomputed` needs its inputs.** It builds nothing, so it requires `--graph`
and one of `--vcfs` or `--graph_alignments` <span class="src">`main.nf:53-56`</span>.

**`-resume`** is Nextflow's own restart and covers the common case of a crashed run; the entry
flags are for runs whose work directory is gone. See [Resuming and skipping work](../guides/skipping-work.md).

---

## Where to read next

- [Stage 1: SV search](../guides/discovery.md) for the callers and the merge.
- [Stage 2: TE annotation](../guides/annotation.md) for RepeatMasker, ULTRA, TSDs and the filters.
- [Stage 3.1: genotyping](../guides/genotyping.md) for PanGenie, Giraffe and GraphAligner.
- [Parameters](../reference/parameters.md) for every flag and its default.
