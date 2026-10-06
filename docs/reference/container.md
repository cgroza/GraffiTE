---
title: Container
description: >-
  The GraffiTE image, the profiles that select it, what is inside it, and where the
  build recipe and the published image are known to differ.
---

# Container

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `25e417a`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## The image

Every process runs in `docker://cgroza/graffite:latest`, pulled from Docker Hub, except
`pav_asm`. All three profiles name the same image. `standard` and `cloud` differ
only in the executor; `cluster` also sets `process.scratch = '$SLURM_TMPDIR'`.
<span class="src">`nextflow.config:7-23`</span>

The image is `linux/amd64` only, 39 layers and about 2.8 GB compressed. Apptainer converts
it to a SIF on first use. On an arm64 host it runs under emulation.

| Profile | Executor | Extra |
|---|---|---|
| `standard` | `local` | |
| `cluster` | `slurm` | `process.scratch = '$SLURM_TMPDIR'` |
| `cloud` | `aws` | `aws` is not a Nextflow executor name; the profile stops at launch. See [Resources and scaling](../guides/resources.md). |

Singularity (or Apptainer) is enabled globally, with automatic mounts. The run options are built
after the params block, so that `--container_tmp` can be read:
<span class="src">`nextflow.config:3-4,176`</span>

```
--contain --bind <--container_tmp, or $(pwd)>:/tmp
```

`--container_tmp` is unset by default, and the bind source is then `$(pwd)`. Pass
`--container_tmp /scratch/you/tmp` and that directory is bound to `/tmp` instead.
<span class="src">`nextflow.config:47`</span> See
[`--container_tmp`](parameters.md#global-switches) and
[Paths and `/tmp`](../getting-started/installation.md#paths-and-tmp).

`--contain` hides the host filesystem from the task except for what is bound in: the task
directory, the input files Nextflow mounts on its own, and `/tmp`, which the bind points at the
task directory or at `--container_tmp`. Two consequences:

- A tool inside a task can only reach files the pipeline staged. Anything a script opens by an
  absolute host path that was not an input is invisible.
- With `--container_tmp` unset, tools that write to `/tmp` write into the task directory, which
  sits on whatever filesystem Nextflow's `work/` is on. Put `work/` on fast, roomy storage, or
  point `--container_tmp` at storage that suits the traffic.

To run it under Docker rather than Apptainer, turn `singularity.enabled` off in a config of
your own and pass `-with-docker cgroza/graffite:latest`.

## What is inside

Built from [`Dockerfile`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/Dockerfile) on top of
`dfam/tetools:latest`, which is Debian 12. The recipe pins few versions; most tools are whatever
the upstream repository, or the base image, held on the day the image was built. We read the last
column out of the published image, `sha256:47836cda` (pushed 2026-09-30), on 2026-10-03.
<span class="src">`Dockerfile:2`</span>

| Tool | Used by | Installed from | Pinned | In `sha256:47836cda` |
|---|---|---|---|---|
| RepeatMasker, RMBlast, HMMER, TRF, RepeatModeler and the rest of the Dfam TETools bundle | `repmask_vcf.sh`, `hervk_ref_state.py` | the `dfam/tetools:latest` base image | no | RepeatMasker 4.2.4, RMBlast 2.17.1+, HMMER 3.4, TRF 4.09, RepeatModeler 2.0.9 <span class="src">`Dockerfile:2`</span> |
| Dfam | nothing: both RepeatMasker calls take `--TE_library` through `-lib` | the base image | no | `dfam40.0.h5` and `dfam40.curated.consensus.0.h5` in `/opt/FamDB-Dfam-4.0/Libraries/famdb` <span class="src">`bin/repmask_vcf.sh:75`, `bin/hervk_ref_state.py:156`</span> |
| Winnowmap, meryl | `map_asm`, `map_longreads` with `--aligner winnowmap` | git, default branch | no | Winnowmap 2.03, meryl r992 <span class="src">`Dockerfile:34-43`</span> |
| htslib (with `bgzip` and `tabix`), samtools, bcftools | nearly every process | git, default branch; bcftools with `--enable-libgsl --enable-perl-filters` | no | htslib 1.24-87-g4be16390, samtools 1.24-47-gde749a60, bcftools 1.24-24-gedf7fd96 <span class="src">`Dockerfile:45-77`</span> |
| SURVIVOR | nothing in v1.1 | git, default branch | no | 1.0.7 <span class="src">`Dockerfile:82-87`</span> |
| minimap2 | `map_asm`, `map_longreads` | git, default branch | no | 2.31-r1302 <span class="src">`Dockerfile:89-94`</span> |
| minigraph | nothing in GraffiTE; panmethyl's `align_minigraph` calls it, and `main.nf` does not import that process | git, default branch | no | 0.21-r606 <span class="src">`Dockerfile:96-101`, `main.nf:12`</span> |
| ULTRA | `repmask_vcf.sh` | git tag | `v1.0.0` | the binary prints no version <span class="src">`Dockerfile:103-111`</span> |
| PanGenie, PanGenie-index | `pangenie_index`, `pangenie` | git, default branch, with cereal `v1.3.2`; the commit is written to `/metadata/pangenie.git.version` | no | commit `d44c4da` <span class="src">`Dockerfile:114-140`</span> |
| pysam, pyparsing, svim-asm, pandas, polars, vcfpy, sniffles, cigar, truvari, pyfaidx, h5py | `svim_asm`, `sniffles_*`, `truvari_merge`, `bin/*.py`, and panmethyl's `annotate_vcf.py` and `merge_csvs.py` | pip, into Debian's Python 3.11 | no | truvari 5.4.0, vcfpy 0.14.2, pysam 0.24.1, pandas 3.0.6, polars 1.44.2, numpy 2.4.6 (pulled in as a dependency), sniffles 2.8.1, svim-asm 1.0.3, pyfaidx 0.9.0.4, pyparsing 3.0.9, cigar 0.1.3, h5py 3.16.0 <span class="src">`Dockerfile:146`</span> |
| R with XML, dplyr, stringr, tidyr, readr, vcfR, optparse | `annotate_vcf.R` | CRAN, into Debian's R | no | R 4.2.2; XML 3.99.0.24, dplyr 1.2.1, stringr 1.6.0, tidyr 1.3.2, readr 2.2.0, vcfR 1.16.0, optparse 1.8.2 <span class="src">`Dockerfile:150`</span> |
| vg | `make_graph`, `graph_align_reads`, `vg_call`, `BED_to_graph` | GitHub release binary | `v1.77.0` | v1.77.0 <span class="src">`Dockerfile:182-183`</span> |
| pypy3 | `subset_gaf.py` (its shebang is `/opt/pypy3/bin/pypy3`) | tarball to `/opt/pypy3` | `7.3.17` | 7.3.17, Python 3.10.14 <span class="src">`Dockerfile:185-190`</span> |
| GraphAligner | `graph_align_reads` with `--graph_method graphaligner` | bioconda, through a throwaway Miniforge | no | 1.0.20 <span class="src">`Dockerfile:194-204`</span> |
| tagtobed, lift_mods, lift_offsets, lift_edges | `bamtags_to_BED` runs `tagtobed` and `lift_epigenome` runs `lift_mods`; the other two belong to panmethyl processes GraffiTE does not import | built from the panmethyl repository with Rust 1.85.0 | `f0aa2c0`, the submodule's commit | not recorded; its `tagtobed` has the panmethyl#3 fix, first in `f0aa2c0` <span class="src">`Dockerfile:206-227`</span> |
| bedtools, ncbi-blast+, pigz, gawk, bc, r-base-core | various | Debian 12 apt | no | bedtools 2.30.0, BLAST 2.12.0+, pigz 2.6, gawk 5.2.1, bc 1.07.1 <span class="src">`Dockerfile:6-31,175`</span> |

## The PAV container

`pav_asm` runs in `library://becklab/pav/pav:latest`, PAV's own image, with 32 CPUs unless
`--cores` says otherwise. GraffiTE's image does not contain PAV. <span class="src">`nextflow.config:339-344`, `module/main.nf:166`</span>

## The recipe and the published image

The published image matches `Dockerfile`: Debian 12, vg `v1.77.0` and the `lift_*` binaries, none
of which `GraffiTE.def` gives. `GraffiTE.def` is the Apptainer recipe that came before it, on
Ubuntu 20.04. We have not built it; read from the recipe, an image built from it would differ from
the published one:

- vg `v1.70.0` rather than `v1.77.0`, and RepeatMasker set up from the TETools sources with the
  Dfam 3.8 root partition <span class="src">`GraffiTE.def:27-164,284`</span>.
- No `lift_mods` and no polars. `lift_epigenome` runs `lift_mods`, and `merge_CSV` runs
  `merge_csvs.py`, which imports polars, so we expect every `--epigenomes` run to stop at one of
  the two <span class="src">`GraffiTE.def:302-310`, `panmethyl/module/main.nf:163,185`</span>.

Two more things follow from the unpinned installs:

- bcftools is built from the default branch. The `--human` filter relies on how bcftools evaluates
  `~` on `Number=.` INFO fields and was checked on bcftools 1.22; the published image has
  1.24-24-gedf7fd96, and a different bcftools may keep different records without an error. See
  [VCF fields](vcf-fields.md).
- `GraffiTE.def` pins numpy to `1.21` and runs `pip3 check`, so a dependency conflict fails its
  build. `Dockerfile` does neither, and the published image has numpy 2.4.6.
  <span class="src">`GraffiTE.def:247-250`, `Dockerfile:146`</span>

To read the versions out of another image:

```bash
apptainer exec graffite.sif bash -c 'vg version | head -1; bcftools --version | head -1; \
  RepeatMasker -v; pip3 list 2>/dev/null | grep -i -E "truvari|vcfpy|pysam|numpy"'
```

## Building it yourself

`Dockerfile` needs Docker with BuildKit, the default builder since Docker Engine 23, because it
uses heredoc `RUN` blocks:

```bash
git clone --recurse-submodules https://github.com/cgroza/GraffiTE.git
cd GraffiTE
docker build -t graffite:local .
```

For Apptainer, convert the result to a SIF:

```bash
apptainer build graffite.sif docker-daemon://graffite:local
```

Then point the profiles at the file instead of the Docker Hub image, in a config passed with `-c`:

```groovy
process.container = '/path/to/graffite.sif'
```

The build pulls `dfam/tetools:latest`, clones several repositories and downloads a Miniforge
installer, so it needs network access and takes a while. A rebuilt image may not carry the same
tool versions as the published one for anything the table above lists as not pinned.

`sudo apptainer build graffite.sif GraffiTE.def` builds the older recipe, with the differences
listed in the previous section.
