---
title: Container
description: >-
  The GraffiTE image, the profiles that select it, what is inside it, and where the
  build recipe and the published image are known to differ.
---

# Container

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## The image

Every process runs in `library://cgroza/collection/graffite:latest`, pulled from the Sylabs
library, except `pav_asm`. All three profiles name the same image; they differ only in the
executor. <span class="src">`nextflow.config:8-25`</span>

| Profile | Executor | Extra |
|---|---|---|
| `standard` | `local` | |
| `cluster` | `slurm` | `process.scratch = '$SLURM_TMPDIR'` |
| `cloud` | `aws` | `aws` is not a Nextflow executor name; the profile stops at launch. See [Resources and scaling](../guides/resources.md). |

Singularity (or Apptainer) is enabled globally, with automatic mounts and these run options:
<span class="src">`nextflow.config:3-5`</span>

```
--contain --bind $(pwd):/tmp
```

`--contain` hides the host filesystem from the task except for what is bound in: the task
directory, the input files Nextflow mounts on its own, and `/tmp`, which the bind points at the
task directory. Two consequences:

- A tool inside a task can only reach files the pipeline staged. Anything a script opens by an
  absolute host path that was not an input is invisible.
- Tools that write to `/tmp` write into the task directory, which sits on whatever filesystem
  Nextflow's `work/` is on. Put `work/` on fast, roomy storage.

The image is also on Docker Hub as `cgroza/graffite`. Nextflow can use it with
`-with-docker cgroza/graffite` once `singularity.enabled` is turned off in your own config.

## What is inside

Built from [`GraffiTE.def`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/GraffiTE.def) on
Ubuntu 20.04. The recipe pins few versions; most tools are whatever the upstream repository's
default branch held on the day the image was built. <span class="src">`GraffiTE.def:1-2`</span>

| Tool | Used by | Installed from | Version |
|---|---|---|---|
| RepeatMasker, RMBlast, HMMER, TRF, RepeatModeler and the rest of the Dfam TETools bundle | `repmask_vcf.sh`, `hervk_ref_state.py` | `TETools/getsrc.sh` at build time | whatever TETools fetched; not recorded <span class="src">`GraffiTE.def:27-165`</span> |
| Dfam library | RepeatMasker | `dfam38_full.0.h5` from the TETools sources | Dfam 3.8 <span class="src">`GraffiTE.def:143-144`</span> |
| Winnowmap, meryl | `map_asm`, `map_longreads` with `--aligner winnowmap` | git, default branch | unpinned <span class="src">`GraffiTE.def:170-175`</span> |
| htslib, samtools, bcftools | nearly every process | git, default branch; bcftools with `--enable-libgsl --enable-perl-filters` | unpinned <span class="src">`GraffiTE.def:177-205`</span> |
| SURVIVOR | nothing in v1.1 | git | unpinned <span class="src">`GraffiTE.def:207-212`</span> |
| minimap2 | `map_asm`, `map_longreads` | git, default branch | unpinned <span class="src">`GraffiTE.def:214-219`</span> |
| ULTRA | `repmask_vcf.sh` | git tag | `v1.0.0` <span class="src">`GraffiTE.def:221-229`</span> |
| PanGenie, PanGenie-index | `pangenie_index`, `pangenie` | git, default branch; the commit is written to `/metadata/pangenie.git.version` in the image | unpinned <span class="src">`GraffiTE.def:233-246`</span> |
| numpy | Python scripts | pip | `1.21` <span class="src">`GraffiTE.def:248`</span> |
| pysam, pyparsing, svim-asm, pandas, vcfpy, sniffles, cigar, truvari, pyfaidx | `svim_asm`, `sniffles_*`, `truvari_merge`, `bin/*.py` | pip | unpinned <span class="src">`GraffiTE.def:250`</span> |
| R with XML, dplyr, stringr, tidyr, readr, vcfR, optparse | `annotate_vcf.R` | CRAN | unpinned <span class="src">`GraffiTE.def:254`</span> |
| vg | `make_graph`, `graph_align_reads`, `vg_call`, `BED_to_graph` | GitHub release binary | `v1.70.0` <span class="src">`GraffiTE.def:285-286`</span> |
| pypy3 | `subset_gaf.py` (its shebang is `/opt/pypy3/bin/pypy3`) | tarball to `/opt/pypy3` | `7.3.17` (Python 3.10) <span class="src">`GraffiTE.def:288-293`</span> |
| GraphAligner | `graph_align_reads` with `--graph_method graphaligner` | bioconda through a throwaway Miniconda | unpinned <span class="src">`GraffiTE.def:295-300`</span> |
| tagtobed | `bamtags_to_BED` | built from the panmethyl repository, default branch | unpinned <span class="src">`GraffiTE.def:303-308`</span> |
| bedtools, ncbi-blast+, tabix, pigz, bc, python3-h5py, r-base-core | various | Ubuntu 20.04 apt | distribution versions <span class="src">`GraffiTE.def:9-25,280`</span> |

## The PAV container

`pav_asm` runs in `library://becklab/pav/pav:latest`, PAV's own image, with 32 CPUs unless
`--cores` says otherwise. GraffiTE's image does not contain PAV. <span class="src">`nextflow.config:321-326`, `module/main.nf:135`</span>

## The recipe and the published image

The recipe in this branch was amended without a rebuild. Before this branch, `GraffiTE.def` did
not install truvari, ULTRA, pyfaidx or pypy3, although the pipeline calls all four and the
published image runs. The additions (ULTRA `v1.0.0`, truvari and pyfaidx from pip, pypy3 `7.3.17`)
are our guess at what the image holds. The amended recipe has not been built, and the versions
in the published image have not been read back. To check them:

```bash
apptainer exec graffite.sif bash -c 'ultra --version; truvari version; \
  python3 -c "import pyfaidx; print(pyfaidx.__version__)"; /opt/pypy3/bin/pypy3 --version'
```

Two more things follow from the unpinned installs:

- bcftools is built from the default branch. The `--human` filter relies on how bcftools evaluates
  `~` on `Number=.` INFO fields and was checked on bcftools 1.22; a different bcftools may keep
  different records without an error. See [VCF fields](vcf-fields.md).
- `pip3 install` without pins can move numpy off `1.21` when truvari asks for a newer one. The
  recipe now runs `pip3 check` after the install so a conflict fails the build instead of the run.
  <span class="src">`GraffiTE.def:248-251`</span>

## Building it yourself

On a Linux host with Apptainer:

```bash
git clone --recurse-submodules https://github.com/cgroza/GraffiTE.git
cd GraffiTE
sudo apptainer build graffite.sif GraffiTE.def
```

Then point the profiles at the file instead of the library, in a config passed with `-c`:

```groovy
process.container = '/path/to/graffite.sif'
```

The build fetches the TETools sources, clones several repositories and downloads a Miniconda
installer, so it needs network access and takes a while. A rebuilt image will not carry the same
tool versions as the published one for anything listed as unpinned above.
