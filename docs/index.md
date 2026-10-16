---
title: GraffiTE
description: >-
  GraffiTE finds transposable element insertion polymorphisms in genome assemblies
  and long-read datasets, annotates them, and genotypes them in read sets using a
  pangenome graph.
---

<h1 class="gt-masthead"><img src="assets/graffite-logo.png" alt="GraffiTE" width="560"></h1>

<p class="gt-tagline">Pangenomic toolbox for the analysis of transposable element insertion polymorphisms.</p>

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `e687712`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](getting-started/v1.0-vs-v1.1.md).

GraffiTE takes genome assemblies or long-read datasets, finds the structural variants that are
transposable element (TE) insertions or deletions relative to a reference genome, annotates what
those elements are, and then genotypes the polymorphisms in as many read sets as you like by
building a pangenome graph with the TEs as bubbles.

It handles both **non-reference** insertions (a TE present in your sample, absent from the
reference) and **reference** insertions (a TE present in the reference, absent from your sample).
You can also use it only as an annotator: hand it a VCF of structural variants and it reports
which of them are TEs.

---

## GraffiTE workflow

A GraffiTE run covers some or all of three stages in series. The third has two parts: genotyping,
and an optional epigenetic mapping. Which stages run depends on the inputs you supply.

```mermaid
%%{init: {"flowchart": {"nodeSpacing": 22, "rankSpacing": 24, "padding": 10}}}%%
flowchart TB
    subgraph S1["Stage 1 · SV search"]
        direction TB
        i_lr(["Long reads<br/>FASTQ or BAM"])
        i_asm(["Genome assemblies"])
        i_svs(["Your SV calls,<br/>one VCF per sample"])
        snf["Map reads,<br/>call SVs<br/><small>Sniffles2</small>"]
        svim["Align assemblies,<br/>call SVs<br/><small>svim-asm</small>"]
        pavc["Call SVs per<br/>haplotype pair<br/><small>PAV</small>"]
        mrg["Merge and collapse<br/>redundant calls<br/><small>truvari collapse</small>"]
        i_lr -->|"<code>--longreads</code><br/><code>--bams</code>"| snf
        i_asm -->|"<code>--assemblies</code>"| svim
        i_asm -->|"<code>--pav</code>"| pavc
        i_svs -->|"<code>--svs</code>"| mrg
        snf --> mrg
        svim --> mrg
        pavc --> mrg
    end
    subgraph S2["Stage 2 · TE annotation"]
        direction TB
        i_vcf(["One merged<br/>SV VCF"])
        rmk["Annotate each SV<br/><small>RepeatMasker, ULTRA</small>"]
        filt["Keep SVs mostly<br/>covered by repeats"]
        hall["Look for hallmarks<br/>of transposition<br/><small>TSD, polyA tail</small>"]
        f_sub[/"pangenome.trusted.vcf, or<br/>pangenome.human.vcf with <code>--human</code>"/]
        i_vcf -.->|"<code>--vcf</code> skips Stage 1"| rmk
        rmk --> filt --> hall --> f_sub
    end
    subgraph S3["Stage 3 · Graph mapping"]
        direction TB
        bld["Build a variation graph<br/><small>each TE polymorphism<br/>is a bubble</small>"]
        i_rd(["Reads to genotype"])
        subgraph S31["3.1 · Genotyping"]
            gt["Genotype each sample<br/><small>vg giraffe or GraphAligner,<br/>then vg call; or PanGenie</small>"]
            gm["Merge the samples<br/>into one VCF"]
            hk["Consolidate HERV-K<br/>loci <small><code>--human</code></small>"]
            gt --> gm -.-> hk
        end
        subgraph S32["3.2 · Epigenetic mapping"]
            lift["Lift methylation<br/>onto the graph<br/><small>panmethyl</small>"]
        end
        bld --> gt
        i_rd -->|"<code>--genotype_with</code>"| gt
        gt -.->|"alignments"| lift
        lift -.-> gm
    end
    out[/"GraffiTE.merged.genotypes.vcf.gz"/]
    S1 -->|"SVs.vcf"| S2
    S2 -->|"pangenome.vcf"| S3
    S3 --> out

    classDef inp fill:#8888881a,stroke:#888888;
    classDef file fill:none,stroke:#888888;
    classDef s1 fill:#f5851f2e,stroke:#f5851f;
    classDef s2 fill:#f2523f2e,stroke:#f2523f;
    classDef s3 fill:#3d63c92e,stroke:#3d63c9;
    classDef s32 fill:#2fb3b32e,stroke:#2fb3b3;
    classDef opt stroke-dasharray:4 3;
    class i_lr,i_asm,i_svs,i_vcf,i_rd inp;
    class out,f_sub file;
    class snf,svim,pavc,mrg s1;
    class rmk,filt,hall s2;
    class bld,gt,gm,hk s3;
    class lift s32;
    class hk,lift opt;
    style S1 fill:#f5851f0d,stroke:#f5851f;
    style S2 fill:#f2523f0d,stroke:#f2523f;
    style S3 fill:#3d63c908,stroke:#3d63c9;
    style S31 fill:#3d63c90d,stroke:#3d63c9;
    style S32 fill:#2fb3b30d,stroke:#2fb3b3,stroke-dasharray:4 3;
```

<p class="gt-legend">Rounded: your input. Slanted: a file GraffiTE writes. Dashed: optional, or a shortcut.
Every stage also reads the reference genome (<code>--reference</code>), and Stage 2 reads a TE
library in RepeatMasker format (<code>--TE_library</code>).</p>

**Stage 1, SV search.** GraffiTE aligns each assembly or read set to the reference, calls
structural variants, and keeps only insertions and deletions. You can mix callers and input types
in one run; GraffiTE merges everything into a single non-redundant SV set.

**Stage 2, TE annotation.** RepeatMasker scans every candidate SV allele against your TE library,
and ULTRA scans it for tandem repeats. GraffiTE discards variants that are not mostly repeat and
annotates the rest with their repeat identity, target site duplications and polyA tails.

**Stage 3.1, genotyping.** GraffiTE builds a pangenome graph with the annotated polymorphisms as
bubbles and genotypes each sample at those in the graph. The default back end,
`--graph_method pangenie`, counts k-mers against a PanGenie index and never aligns; `giraffe`
and `graphaligner` map the reads onto the graph and call genotypes with `vg call`.

**Stage 3.2, epigenetic mapping, optional.** With `--epigenomes`, on the `giraffe`,
`graphaligner` and `precomputed` graph methods (every method but `pangenie`), the `panmethyl`
submodule lifts the base modifications carried by the genotyping BAMs onto the graph and
summarizes them per polymorphism. See [Epigenetic mapping](guides/methylation.md).

Stage 3 is optional (`--genotype false`). Stages 1 and 2 can each be skipped by supplying their
outputs directly; see [Resuming and skipping work](guides/skipping-work.md).

---

## Reference and non-reference insertions

GraffiTE reports variants **relative to the reference genome**. Read the direction off the
allele lengths: ALT longer than REF is a
non-reference insertion, REF longer than ALT is a reference insertion. `SVLEN` is positive for the
first and negative for the second. Do not use `SVTYPE` for this. On any run with two or more
caller VCFs, `truvari_merge` strips every upstream INFO field before the collapse and adds back
only `SVLEN`, so `pangenome.vcf` has no `SVTYPE`. It survives only on the single-input paths,
`--vcf` and a single caller VCF.

<figure>
--8<-- "assets/ref-vs-nonref-insertion.svg"
<figcaption>
ALT longer than REF means the element is <em>absent</em> from the reference and
<em>present</em> in the sample. REF longer than ALT means the element is <em>present</em> in
the reference and <em>absent</em> from the sample. The <code>SVTYPE</code> values shown appear
only on the single-input paths. In both cases GraffiTE annotates the element itself.
</figcaption>
</figure>

**An ALT allele does not mean "TE present"**. For a reference insertion the relationship is
inverted. `vcf_to_pa_tsv.py` writes the presence-absence TSVs and does that conversion for you:
a `1` always means the TE is in that sample, whichever way the record points. It reads the
direction from `INFO/SVTYPE` where there is one and from the allele lengths where there is not,
which is what a run with two or more caller VCFs produces. It writes `NA` where neither gives a
direction: on a symbolic or multi-allelic ALT, or when REF and ALT are the same length. See
[Output files](reference/outputs.md).

---

## Glossary

Where the codebase uses a term loosely, this table states the meaning that applies here.

| Term | Meaning |
|---|---|
| **MEI** | Mobile element insertion; on these pages, one that is present in some haplotypes and absent from others. The v1.0 README calls them pMEs, polymorphic mobile elements. We use MEI, the more common term. |
| **Non-reference insertion** | TE present in the sample, absent from the reference. ALT longer than REF, `SVLEN` positive. |
| **Reference insertion** | TE present in the reference, absent from the sample. REF longer than ALT, `SVLEN` negative. |
| **Hit** | One RepeatMasker match after fragment grouping, a single element. Counted by the `n_hits` INFO field. |
| **Fragment** | One raw line of RepeatMasker output. Several fragments may be grouped into one hit. Counted by `fragmts`. |
| **Repeat span** | The fraction of a variant's sequence covered by the non-redundant union of RepeatMasker TE hits and ULTRA tandem repeats. The `total_repeat_span` INFO field; the main quality filter. |
| **Trusted subset** | A conservative subset of `pangenome.vcf`: single-hit, at least 250 bp, not dominated by tandem repeat, and polyA-supported if it is a non-LTR retrotransposon. Written to `pangenome.trusted.vcf`. |
| **Human MEI subset** | With `--human`, a subset filtered to recent human mobile element subfamilies (AluY, L1HS, SVA_D/E/F, HML-2). Written to `pangenome.human.vcf`, **instead of** the trusted subset. |
| **TSD** | Target site duplication. A short direct repeat flanking a genuine mobile element insertion, created by the integration mechanism. |
| **Locus** (HERV-K) | With `--human`, the set of records that describe the same HML-2 element. Two records share a locus when they touch the same reference HML-2 element, or when their footprints lie within `--hervk_locus_window` (`1200` bp) of each other. Named by the `HERVK_LOCUS` INFO field. |
| **Consolidated VCF** | With `--human`, a VCF in which each HERV-K locus is one record with per-allele states (`pangenome.human.consolidated.vcf`, and `GraffiTE.merged.genotypes.human.vcf.gz` after genotyping). |
| **Precomputed** | The `--graph_method precomputed` mode: no graph is built and no reads are mapped; the run reuses a `--graph` directory with `--vcfs` or `--graph_alignments` from an earlier run. |
| **Stage 1 / 2 / 3** | SV search / TE annotation / graph mapping (3.1 genotyping, 3.2 epigenetic mapping), as above. |

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

- :material-flask-outline: **[Epigenetic mapping](guides/methylation.md)**

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
