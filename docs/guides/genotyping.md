---
title: Stage C: genotyping
description: >-
  How GraffiTE builds a pangenome graph from the annotated polymorphisms, maps read sets onto
  it with PanGenie, Giraffe or GraphAligner, and merges the per-sample genotypes.
---

# Stage C: genotyping

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

Stage C takes `pangenome.vcf` from [Stage B](annotation.md) and the read sets listed in
`--genotype_with`, builds a graph in which every polymorphism is a bubble, and genotypes each
sample at each bubble. It runs by default; `--genotype false` stops after Stage B
<span class="src">`main.nf:188`</span>.

---

## What a bubble is

<figure>
--8<-- "assets/graph-bubble.svg"
<figcaption>
Each polymorphism becomes two paths between the same flanks. For an <code>INS</code> record the
reference path skips the element and the ALT path runs through it. For a <code>DEL</code> record
it is the other way round. Each read counts for whichever path it traverses, so reads through
the TE node always come from haplotypes that carry the element, whichever allele that is.
</figcaption>
</figure>

The consequence for reading the output is the one stated on the [home page](../index.md): the
ALT allele means *TE present* for `INS` records and *TE absent* for `DEL` records.

---

## Read sets

`--genotype_with` is a samplesheet with `sample`, `path` and `type` columns
(see [Samplesheets](../reference/samplesheets.md)). The `type` column selects a Giraffe
parameter preset <span class="src">`main.nf:189-205`</span>:

| `type` | preset | used by |
|---|---|---|
| `pb`, `hifi` | `hifi` | `vg giraffe --parameter-preset hifi` |
| `ont` | `r10` | `vg giraffe --parameter-preset r10` |
| anything else, including an empty column | `default` | `vg giraffe --parameter-preset default -i` |

Two details follow from that table:

- With the `default` preset, `vg giraffe` is given `-i`, so the FASTQ is read as
  **interleaved paired-end** <span class="src">`module/main.nf:773-776, 781`</span>. Short-read
  samples must be supplied as a single interleaved file, not as two mate files.
- PanGenie ignores the preset: it counts k-mers and never aligns
  <span class="src">`module/main.nf:699`</span>.

A `path` ending in `.bam` goes through `bam_to_fastq` first: alignment tags are stripped, the
file is name-sorted and converted back to FASTQ with `samtools fastq`
<span class="src">`main.nf:206-211`, `module/main.nf:744-760`</span>. The alignments in the BAM
are not used; only the reads are. Methylation tags in such a BAM are read separately, see
[Methylation](methylation.md).

---

## Choosing a graph method

| `--graph_method` | Graph | Reads are | Genotyper | When |
|---|---|---|---|---|
| `pangenie` (default) | `PanGenie-index` on a merged VCF | counted as k-mers | PanGenie | k-mer counting, no alignment step |
| `giraffe` | `vg autoindex` (GBZ) | aligned with `vg giraffe` | `vg call` | short or long reads; the back end the HERV-K consolidation is validated on |
| `graphaligner` | `vg construct` (GFA) | aligned with `GraphAligner` | `vg call` | long reads |
| `precomputed` | supplied with `--graph` | supplied with `--graph_alignments`, or skipped with `--vcfs` | `vg call`, or none | re-genotyping an existing graph |

Source: <span class="src">`main.nf:214-275`</span>, <span class="src">`nextflow.config:35`</span>.

Anything else stops the run:

```text
Unsupported --graph_method. --graph_method must be pangenie, giraffe, graphaligner or precomputed.
```

`precomputed` builds nothing itself. Without `--graph` and one of `--vcfs` or
`--graph_alignments` the run stops before any process starts
<span class="src">`main.nf:59-64`</span>:

```text
--graph_method precomputed builds nothing itself: it needs --graph (an index directory holding
index.gfa and index.pb, as make_graph writes) and either --vcfs (per-sample vg call VCFs) or
--graph_alignments (per-sample gaf,pack).
```

See [Resuming and skipping work](skipping-work.md) for what each of those inputs looks like.

---

## PanGenie

Two processes <span class="src">`module/main.nf:667-706`</span>.

**`pangenie_index`**, once per run:

1. `bcftools view -G` drops the genotype columns from `pangenome.vcf`; an `awk` step adds a single
   pseudo-sample `ref` carrying `1|0` at every record, the shape PanGenie's graph builder expects.
2. Records at the same position are joined into multi-allelic sites (`bcftools norm -m+`).
3. `merge_vcfs.py merge -ploidy 2` (PanGenie's own helper, shipped in `bin/`) turns the file into
   the multi-sample graph VCF PanGenie reads, with an `INFO/ID` field naming the GraffiTE variant
   behind each ALT allele.
4. `PanGenie-index` builds the k-mer index.

**`pangenie`**, once per sample: `PanGenie -s <sample> -i <(zcat -f reads) -f pangenie_index`,
then `bcftools norm -f ref -m-` splits the multi-allelic genotypes back into one record per
GraffiTE variant, so that they match `pangenome.vcf` one for one when the merge transfers INFO.
The result is published as `4_Genotyping/<sample>_genotyping.vcf.gz` with its index.

**Resources:** `--pangenie_threads`, `--pangenie_memory`, `--pangenie_time` for both processes
<span class="src">`nextflow.config:235-244`</span>.

---

## Giraffe and GraphAligner

Three processes, plus `bam_to_fastq` when needed.

### `make_graph`

<span class="src">`module/main.nf:708-742`</span>. `bcftools +setGT -- -t a -n u` unphases every genotype
in `pangenome.vcf`, then:

| method | commands | `index/` holds |
|---|---|---|
| `giraffe` | `vg autoindex -w sr-giraffe -w lr-giraffe`, `vg convert --vg-algorithm -f`, `vg snarls` | `index.giraffe.gbz` and its companions, `index.gfa`, `index.pb` |
| `graphaligner` | `vg construct -a -m 1024`, `vg convert --vg-algorithm -f`, `vg snarls` | `index.vg`, `index.gfa`, `index.pb` |

`index.gfa` is the graph as GFA, `index.pb` the snarl (bubble) decomposition `vg call` needs.
The directory is published as `GraffiTE_graph/index/` and is what `--graph` takes on a later
run. Skipped entirely when `--graph` is given <span class="src">`main.nf:221-225`</span>.

**Resources:** `--make_graph_threads`, `--make_graph_memory` (default `40G`), `--make_graph_time`
(default `6h`) <span class="src">`nextflow.config:245-249`</span>.

### `graph_align_reads`

<span class="src">`module/main.nf:762-798`</span>, once per sample:

| method | aligner | then |
|---|---|---|
| `giraffe` | `vg giraffe --parameter-preset <preset> -o gam --index-basename index/index [-i] -f reads` | `vg pack -Q <min_mapq>` and `vg convert -G` to GAF |
| `graphaligner` | `GraphAligner -x vg -g index/index.gfa -f reads -a sample.gam` | same |

`vg pack` builds the per-node coverage `vg call` reads, dropping alignments below
`--min_mapq` (default `0`) <span class="src">`nextflow.config:120`</span>. The GAF is passed
through `subset_gaf.py`, which keeps the twelve standard columns plus the `cs`/`cg` difference
string, then sorted by read name and gzipped. The GAM is deleted. Published as
`GraffiTE_alignments/<sample>.gaf.gz` and `<sample>.pack`, which is the pair
`--graph_alignments` takes on a later run.

This process runs with `errorStrategy = 'finish'`: a failed sample lets the other samples
complete before the run stops <span class="src">`nextflow.config:255-260`</span>.

**Resources:** `--graph_align_threads`, `--graph_align_memory`, `--graph_align_time` (default
`12h`); `bam_to_fastq` uses the same three.

### `vg_call`

<span class="src">`module/main.nf:800-816`</span>, once per sample:

```bash
vg call -a -A --threads N -R chrX:1,chrY:1 -m 2,4 -r index/index.pb -s <sample> -k <sample>.pack index/<graph> \
  | bcftools norm -m- | bcftools sort -Oz -o <sample>.vcf.gz
```

- `<graph>` is `index.giraffe.gbz` for `giraffe` and `index.gfa` otherwise, `precomputed`
  included.
- `-a` and `-A` request a call at every snarl, reference calls and nested snarls included, so that
  every sample has a genotype at every bubble and the merge lines up. Check `vg call --help` for
  the exact wording in your vg version; the container ships vg 1.70.0.
- `-m 2,4` is `--min_support`, passed to vg as its minimum allele and site support
  <span class="src">`nextflow.config:121`</span>.
- `-R chrX:1,chrY:1` sets ploidy 1 on the contigs named **exactly** `chrX` and `chrY`. Every other
  contig is called diploid, so a reference that names them `X` and `Y`, or `NC_000023.11`, gets
  diploid calls on the sex chromosomes. `vg call` has no other ploidy input.
- `bcftools norm -m-` splits multi-allelic calls into one record per ALT.

**Resources:** `--vg_call_threads`, `--vg_call_memory`, `--vg_call_time` (default `2h`)
<span class="src">`nextflow.config:261-265`</span>.

---

## The merged genotypes

`merge_VCFs` <span class="src">`module/main.nf:818-844`</span> takes every per-sample VCF and:

1. `bcftools merge -m none` joins them into one multi-sample VCF without creating multi-allelic
   records.
2. `bcftools annotate -a pangenome.vcf -c CHROM,POS,ID,REF,ALT,INFO` copies every INFO field of
   the matching `pangenome.vcf` record, so the repeat annotation, TSD and polyA fields travel
   with the genotypes. A genotyped record matches only when `CHROM`, `POS`, `ID`, `REF` and `ALT`
   all agree.
3. The `##GraffiTE_version` header line is added.

Published as `4_Genotyping/GraffiTE.merged.genotypes.vcf.gz`. No `.tbi` is published for it: the
index written earlier in the script belongs to the pre-annotation file, which is deleted. Run
`tabix -p vcf` on it yourself.

**Resources:** `--merge_vcf_memory` (default `10G`), `--merge_vcf_time` (default `1h`)
<span class="src">`nextflow.config:266-270`</span>.

!!! note "No presence-absence TSV for the genotyped calls"
    The presence-absence TSVs in `3_TSD_search/` are built from the discovery genotypes in
    `pangenome.vcf`. Nothing converts `GraffiTE.merged.genotypes.vcf.gz` the same way; apply the
    `INS`/`DEL` polarity rule above yourself, or start from the TSV's column schema in
    [Output files](../reference/outputs.md).

With `--human`, the human subset of this file goes through one more step,
[`hervk_reconcile`](human-mei.md), which writes `GraffiTE.merged.genotypes.human.vcf.gz`.

---

## Reading a genotype

| Record | `GT` | The TE is |
|---|---|---|
| `SVTYPE=INS` | `0/0` | absent from this sample |
| `SVTYPE=INS` | `0/1`, `1/1` | present on one or both haplotypes |
| `SVTYPE=DEL` | `0/0` | present (the reference carries it) |
| `SVTYPE=DEL` | `0/1`, `1/1` | absent from one or both haplotypes |

FORMAT fields beyond `GT` are whatever `vg call` or PanGenie wrote; see
[VCF fields](../reference/vcf-fields.md).

---

## Next

- [Human MEIs](human-mei.md): the HERV-K consolidation that follows genotyping under `--human`
- [Methylation](methylation.md): lifting base-modification calls onto the graph
- [Resuming and skipping work](skipping-work.md): `--graph`, `--graph_alignments`, `--vcfs`
- [Output files](../reference/outputs.md)
