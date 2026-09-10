---
title: Methylation
description: >-
  Lifting base-modification calls from long-read BAMs onto the pangenome graph with the
  panmethyl submodule, and reading methylation levels per allele in the genotyped VCF.
---

# Methylation

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

!!! warning "Documented from the wiring, not from a run"
    This page describes what `main.nf` and the panmethyl module do when `--epigenomes` is set.
    No methylation run was executed while writing it. Process behaviour is read from
    `panmethyl/module/main.nf` at submodule commit `bd1c383`.

With `--epigenomes`, GraffiTE reads the base-modification tags (`MM`/`ML`) of the BAMs it is
genotyping, lifts each modified base onto the graph node it maps to, and reports per-allele
methylation levels in the genotyped VCF. The processes come from
[panmethyl](https://github.com/cgroza/panmethyl), a git submodule that must be present for
`main.nf` to parse at all; see [Installation](../getting-started/installation.md).

---

## Requirements

| Requirement | Why | Source |
|---|---|---|
| `--graph_method giraffe`, `graphaligner` or `precomputed` | the methylation branch sits inside the vg genotyping block; PanGenie has no graph alignments to lift onto | <span class="src">`main.nf:217, 245`</span> |
| `--epigenomes` | switches the branch on | <span class="src">`nextflow.config:30`</span> |
| BAM entries in `--genotype_with` | modifications are read from BAM tags; a FASTQ sample is genotyped but gets no methylation | <span class="src">`main.nf:255`</span> |
| `MM`/`ML` tags in those BAMs | what `tagtobed` extracts | <span class="src">`panmethyl/module/main.nf:146`</span> |

The reads of a BAM sample are also extracted to FASTQ and aligned to the graph as for any other
sample ([Stage C](genotyping.md)). The alignment inside the BAM is not used; the graph alignment
(`<sample>.gaf.gz`) is what the modifications are lifted through.

---

## How it works

<span class="src">`main.nf:245-268`</span>

1. **`index_graph`** reads `index.gfa` and lists every position of `--motif` (default `CG`) on
   every node, writing `node_sizes.csv`, `nodes_list.csv` and `index.csv.gz` to `out/index/`
   <span class="src">`panmethyl/module/main.nf:67-83`</span>.
2. **`bamtags_to_BED`** runs `tagtobed` on each BAM, extracting the modification calls for
   `--code` (default `C+m`, 5-methylcytosine on the forward strand) into `<sample>.mods.gz`
   <span class="src">`panmethyl/module/main.nf:136-148`</span>.
3. **`lift_epigenome`** joins the modification calls with the sample's graph alignment by read
   name and projects each modified base onto graph node coordinates with `lift_mods`, writing
   `out/lifted/<sample>.csv.gz` <span class="src">`panmethyl/module/main.nf:150-166`</span>.
4. **`merge_CSV`** aggregates the per-read calls into per-node levels at the indexed motif
   positions, writing `out/levels/<sample>.csv.gz`
   <span class="src">`panmethyl/module/main.nf:168-187`</span>.
5. **`annotate_VCF`** walks each allele's path (`INFO/AT`) in the sample's `vg call` VCF, sums
   the levels on the nodes of that path, and writes `out/annotation/<sample>.mods.vcf.gz` plus a
   long-format `<sample>.mods.tsv` with one row per allele
   <span class="src">`panmethyl/module/main.nf:1-15`, `panmethyl/module/resources/usr/bin/annotate_vcf.py:44-46`</span>.

The annotated per-sample VCFs replace the plain `vg call` VCFs going into `merge_VCFs`, so the
merged genotypes carry the methylation FORMAT fields <span class="src">`main.nf:262`</span>.

### The FORMAT fields

Added by `annotate_vcf.py` <span class="src">`panmethyl/module/resources/usr/bin/annotate_vcf.py:44-46`</span>, one value per allele of the genotype:

| Field | Type | Meaning |
|---|---|---|
| `PMN` | Float, `Number=.` | number of modified nucleotides with coverage on the allele's path |
| `PML` | Float, `Number=.` | average modification level across the path |
| `PMD` | Float, `Number=.` | average depth across the modified nucleotides |

`PML` per allele is what lets a TE insertion's methylation be compared against the empty site,
or one haplotype's copy against the other's.

---

## Inputs

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--epigenomes` | `false` | run the branch | <span class="src">`nextflow.config:30`</span> |
| `--motif` | `"CG"` | motif indexed on the graph | <span class="src">`nextflow.config:165`</span> |
| `--code` | `"C+m"` | modification code passed to `tagtobed` (`-T C -B C+m`) | <span class="src">`nextflow.config:164`, `panmethyl/module/main.nf:146`</span> |
| `--lifted` | `false` | samplesheet (`sample,path`) of already-lifted `<sample>.csv.gz` files; skips steps 2 and 3 | <span class="src">`main.nf:250-253`</span> |
| `--bed` | `false` | a BED of regions to project onto the graph and annotate with methylation | <span class="src">`main.nf:264-268`</span> |

With `--bed`, three more processes run: `BED_to_graph` projects the regions with `vg annotate`,
`annotate_BED` sums levels over each region per sample into `out/annotation/<sample>.bed`, and
`merge_BED` joins the samples into `out/annotation/merged_epiannoation.bed` (the filename is
spelled that way in the module) <span class="src">`panmethyl/module/main.nf:17-65`</span>.

---

## Outputs

All relative to `--out`, and beside GraffiTE's numbered directories:

| Path | Written by |
|---|---|
| `index/node_sizes.csv`, `nodes_list.csv`, `index.csv.gz` | `index_graph` |
| `lifted/<sample>.csv.gz` | `lift_epigenome` |
| `levels/<sample>.csv.gz` | `merge_CSV` |
| `annotation/<sample>.mods.vcf.gz`, `<sample>.mods.tsv` | `annotate_VCF` |
| `annotation/<sample>.bed`, `merged_epiannoation.bed` | `annotate_BED`, `merge_BED` (with `--bed`) |
| `4_Genotyping/GraffiTE.merged.genotypes.vcf.gz` | `merge_VCFs`, now carrying `PMN`, `PML`, `PMD` |

---

## Resources

The panmethyl processes have fixed allocations in `nextflow.config` rather than parameters
<span class="src">`nextflow.config:281-320`</span>:

| Process | CPUs | Memory | Time |
|---|---|---|---|
| `bamtags_to_BED` | 2 | 50 GB | 6 h |
| `lift_epigenome`, `merge_CSV` | 1 | 60 GB | 6 h |
| `index_graph`, `annotate_VCF`, `annotate_BED`, `BED_to_graph`, `merge_BED` | 1 | 40 GB | 6 h |

Change them with a `-c` config file that overrides the `withName` blocks.

---

## Limitations

- PanGenie runs get no methylation: the branch needs graph alignments.
- One modification code per run.
- panmethyl's own aligners (`align_giraffe`, `align_graphaligner`, `align_minigraph`) and its
  `lift_nucleotides` process are not called by GraffiTE; the alignment comes from
  `graph_align_reads`.
