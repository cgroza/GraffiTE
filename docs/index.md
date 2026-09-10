---
title: GraffiTE
description: >-
  GraffiTE finds transposable element insertion polymorphisms in genome assemblies
  and long-read datasets, annotates them, and genotypes them in read sets using a
  pangenome graph.
---

# GraffiTE

**Pangenomic toolbox for the analysis of transposable element insertion polymorphisms.**

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](getting-started/v1.0-vs-v1.1.md).

GraffiTE takes genome assemblies or long-read datasets, finds the structural variants that are
transposable element (TE) insertions or deletions relative to a reference genome, annotates what
those elements are, and then genotypes each polymorphism in as many read sets as you like by
building a pangenome graph in which every TE is a bubble.

It handles both **non-reference** insertions (a TE present in your sample, absent from the
reference) and **reference** insertions (a TE present in the reference, absent from your sample).
It can also be used purely as an annotator: hand it a VCF of structural variants and it will tell
you which of them are TEs.

---

## The three stages

Every GraffiTE run is some subset of three stages in series, plus an optional fourth. Which stages
run depends on which inputs you supply.

```mermaid
flowchart LR
    subgraph A["Stage A · Discovery"]
        direction TB
        A1["Assemblies<br/><code>--assemblies</code> · <code>--pav</code>"]
        A2["Long reads<br/><code>--longreads</code> · <code>--bams</code>"]
        A3["Existing SV calls<br/><code>--svs</code> · <code>--vcf</code>"]
        A1 --> AM["Merge<br/><small>truvari collapse</small>"]
        A2 --> AM
        A3 --> AM
    end

    subgraph B["Stage B · Annotation"]
        direction TB
        B1["RepeatMasker + ULTRA<br/><small>what repeat is it?</small>"]
        B2["Repeat-span filter<br/><small>is it mostly repeat?</small>"]
        B3["TSD + polyA<br/><small>does it look like a real MEI?</small>"]
        B1 --> B2 --> B3
    end

    subgraph C["Stage C · Genotyping"]
        direction TB
        C1["Build pangenome graph"]
        C2["Map reads to graph"]
        C3["Call genotypes"]
        C1 --> C2 --> C3
    end

    AM --> B1
    B3 --> C1

    A ~~~ B ~~~ C
```

**Stage A, discovery.** Each assembly or read set is aligned to the reference, structural
variants are called, and only insertions and deletions are kept. Multiple callers and multiple
input types can be mixed in one run; everything is merged into a single non-redundant SV set.

**Stage B, annotation.** Every candidate SV allele is scanned against your TE library with
RepeatMasker and against itself with ULTRA (tandem repeats). Variants that are not mostly repeat
are discarded. Survivors are annotated with their repeat identity, target site duplications and
polyA tails.

**Stage C, genotyping.** The annotated polymorphisms are induced into a pangenome graph as
bubbles, reads are mapped onto it, and each sample is genotyped at every polymorphism.

**Methylation, optional.** With `--epigenomes`, and only on the `giraffe` or `graphaligner`
graph methods, the base modifications carried by the genotyping BAMs are lifted onto the graph
and summarised per polymorphism. This step lives in the `panmethyl` submodule. See
[Methylation](guides/methylation.md).

Stage C is optional (`--genotype false`). Stages A and B can each be skipped by supplying their
outputs directly; see [Resuming and skipping work](guides/skipping-work.md).

---

## Reference and non-reference insertions

GraffiTE reports variants **relative to the reference genome**, so the VCF `SVTYPE` and the
presence of the TE point in opposite directions half the time:

<figure>
--8<-- "assets/ref-vs-nonref-insertion.svg"
<figcaption>
An <code>INS</code> record means the element is <em>absent</em> from the reference and
<em>present</em> in the sample. A <code>DEL</code> record means the element is <em>present</em> in
the reference and <em>absent</em> from the sample. In both cases the element itself is what
GraffiTE annotated.
</figcaption>
</figure>

So **an ALT allele does not mean "TE present"**. For `DEL` records the relationship is inverted.
The presence-absence TSVs published alongside each VCF do this conversion for you: a `1` always
means *the TE is present in this sample*, whichever way the VCF record points. See
[Output files](reference/outputs.md).

---

## Glossary

Where the codebase uses a term loosely, this table states the meaning that applies here.

| Term | Meaning |
|---|---|
| **pME** | Polymorphic mobile element. A mobile element insertion that is present in some haplotypes and absent in others. |
| **Non-reference insertion** | TE present in the sample, absent from the reference. Appears as `SVTYPE=INS`. |
| **Reference insertion** | TE present in the reference, absent from the sample. Appears as `SVTYPE=DEL`. |
| **Hit** | One RepeatMasker match after fragment grouping, a single element. Counted by the `n_hits` INFO field. |
| **Fragment** | One raw line of RepeatMasker output. Several fragments may be grouped into one hit. Counted by `fragmts`. |
| **Repeat span** | The fraction of a variant's sequence covered by the non-redundant union of RepeatMasker TE hits and ULTRA tandem repeats. The `total_repeat_span` INFO field; the main quality filter. |
| **Trusted subset** | A conservative subset of `pangenome.vcf`: single-hit, long enough, not dominated by tandem repeat, and polyA-supported if it is a non-LTR element. Written to `pangenome.trusted.vcf`. |
| **Human pME subset** | With `--human`, a subset filtered to recent human mobile element subfamilies (AluY, L1HS, SVA_D/E/F, HML-2). Written to `pangenome.human.vcf`, **instead of** the trusted subset. |
| **TSD** | Target site duplication. A short direct repeat flanking a genuine mobile element insertion, created by the integration mechanism. |
| **Locus** (HERV-K) | With `--human`, the set of records that describe the same HML-2 element, grouped by overlap. Named by the `HERVK_LOCUS` INFO field. |
| **Consolidated VCF** | With `--human`, a VCF in which each HERV-K locus is one record with per-allele states (`pangenome.human.consolidated.vcf`, and `GraffiTE.merged.genotypes.human.vcf.gz` after genotyping). |
| **Precomputed** | The `--graph_method precomputed` mode: no graph is built and no reads are mapped; the run reuses a `--graph` directory with `--vcfs` or `--graph_alignments` from an earlier run. |
| **Stage A / B / C** | Discovery / annotation / genotyping, as above. |

---

## Where to go next

<div class="grid cards" markdown>

- :material-download: **[Installation](getting-started/installation.md)**

    Prerequisites, cloning with submodules, and the container image.

- :material-rocket-launch: **[Quickstart](getting-started/quickstart.md)**

    Run the bundled human chromosome-22 test dataset end to end.

- :material-sign-direction: **[Choosing your inputs](getting-started/choosing-your-inputs.md)**

    A decision tree from "what data do I have" to "which flags do I pass".

- :material-tune: **[Parameters](reference/parameters.md)**

    Every parameter, its default, and what it controls.

- :material-file-tree: **[Output files](reference/outputs.md)**

    What GraffiTE writes and where.

- :material-dna: **[Human MEIs](guides/human-mei.md)**

    The `--human` filter and the HERV-K classifier.

- :material-flask-outline: **[Methylation](guides/methylation.md)**

    Lift base modifications from the genotyping BAMs onto the graph with `--epigenomes`.

- :material-compare: **[v1.0 vs v1.1](getting-started/v1.0-vs-v1.1.md)**

    What changed since the paper, flag by flag and field by field.

</div>

---

## Citing GraffiTE

> Groza, C., Chen, X., Wheeler, T.J. et al. A unified framework to analyze transposable element
> insertion polymorphisms using graph genomes. *Nature Communications* **15**, 8915 (2024).
> [doi:10.1038/s41467-024-53294-2](https://doi.org/10.1038/s41467-024-53294-2)

GraffiTE was developed by **Cristian Groza** and **Clément Goubert** in
[Guillaume Bourque's group](https://computationalgenomics.ca/BourqueLab/) at the
[McGill Genome Centre](https://www.mcgillgenomecentre.ca/), Montréal, Canada. It builds on the
concept described in [Groza et al., 2022](https://link.springer.com/protocol/10.1007/978-1-0716-2883-6_5).

Bugs, comments and suggestions are welcome in the
[issue tracker](https://github.com/cgroza/GraffiTE/issues).
