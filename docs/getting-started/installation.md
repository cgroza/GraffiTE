---
title: Installation
description: >-
  What GraffiTE needs on the machine, how to get the pipeline and its container
  image, and how to check the install before a real run.
---

# Installation

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1.0.0` at commit `68ddc42`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](v1.0-vs-v1.1.md).

GraffiTE is a Nextflow pipeline. Every tool it calls lives in a container image, so the host
needs Nextflow, a Java runtime for it, and a container engine. It runs on any Linux machine or
cluster; the `cluster` profile targets SLURM.

---

## Prerequisites

| Requirement | Notes |
|---|---|
| Java 17 or newer | Nextflow's own requirement. |
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | See the note below on the parser. |
| [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or Singularity | The container engine the config enables by default. Docker works with a config change. |

!!! note "Nextflow 26"
    Nextflow 26.04 parses scripts with its strict syntax by default. `main.nf` compiles under it
    since commit `4605630` (2026-09-10, merged into `v1.1dev` on 2026-09-15 in PR #101); previewed
    here with 26.04.6 under both parsers. Checkouts older than that commit need
    `export NXF_SYNTAX_PARSER=v1`.

!!! note "Apptainer through Conda"
    Users have reported problems with Apptainer installed through Conda. Install it from the
    Apptainer packages instead.

---

## Get GraffiTE

Two ways, depending on whether you want a local copy of the code.

=== "Let Nextflow fetch it"

    Nextflow downloads and caches the pipeline the first time you run it, and initializes the
    `panmethyl` submodule for you:

    ```bash
    nextflow pull cgroza/GraffiTE -r v1.1.0.0
    ```

    The cached copy lives under `~/.nextflow/assets/cgroza/GraffiTE/`. Later runs use
    `nextflow run cgroza/GraffiTE -r v1.1.0.0 ...`.

=== "Clone the repository"

    ```bash
    git clone --branch v1.1.0.0 --recurse-submodules https://github.com/cgroza/GraffiTE.git
    cd GraffiTE
    ```

    `main.nf` includes the `panmethyl` submodule unconditionally
    <span class="src">`main.nf:12`</span>, and a clone made without `--recurse-submodules` leaves
    that directory empty, so Nextflow stops at compile time with:

    ```
    Error main.nf:12:1: Invalid include source: '<dir>/panmethyl/module.nf'
    ```

    Run `git submodule update --init` in the clone and start again. `nextflow pull` fetches the
    submodule itself.

---

## Get the container image

The config points every process but `pav_asm` at `docker://cgroza/graffite:1.1.0.0`
<span class="src">`nextflow.config:10,15`</span>, and Nextflow pulls it on first use. Each release
has its own image tag; `latest` follows the newest release.
Apptainer converts the Docker image to a SIF itself, so you have nothing else to register.

The image is `linux/amd64` only. On an arm64 machine it runs under emulation, slowly.

If the compute nodes have no internet access, pull the image once on a node that does and point
the run at the file:

```bash
apptainer pull graffite_1.1.0.0.sif docker://cgroza/graffite:1.1.0.0
nextflow run cgroza/GraffiTE -r v1.1.0.0 -with-singularity /abs/path/graffite_1.1.0.0.sif ...
```

`-with-singularity` overrides the image path in `nextflow.config`.

The PAV entry point (`--pav`) uses a second image, `library://becklab/pav/pav:latest`
<span class="src">`nextflow.config:345`</span>, pulled the same way.
[Container contents](../reference/container.md) lists what the image contains and which recipe it
is built from.

### Paths and `/tmp`

Nextflow runs every container with `--contain --bind <dir>:/tmp`
<span class="src">`nextflow.config:181`</span>. Two consequences:

- The host filesystem is hidden from the container except for what Nextflow mounts: the work
  directory and the files it stages (`singularity.autoMounts = true`
  <span class="src">`nextflow.config:4`</span>). Give input files, including the paths inside
  samplesheets, as absolute paths.
- `/tmp` inside the container is the directory the task runs in, unless you say otherwise: its
  directory under `work/`, or, under `-profile cluster`, its scratch directory in `$SLURM_TMPDIR`
  on the node. Where that filesystem is small, slow, or not writable from inside the container,
  point `--container_tmp` at one that is:

  ```bash
  nextflow run cgroza/GraffiTE -r v1.1.0.0 --container_tmp /scratch/you/tmp ...
  ```

  The directory must exist before the run starts. On SLURM sites that set a per-job scratch
  variable, `--container_tmp $SLURM_TMPDIR` works too.

!!! warning "Do not edit the cached config to do this"
    Before v1.1 the README's fix was to edit `singularity.runOptions` in
    `~/.nextflow/assets/cgroza/GraffiTE/nextflow.config`. That works once and then makes
    `nextflow pull` and `-latest` fail with `contains uncommitted changes`
    <span class="src">[#93](https://github.com/cgroza/GraffiTE/issues/93)</span>. Use
    `--container_tmp`, or `-c` with a config file of your own.

A container that cannot write `/tmp` does not always fail. The run can finish having skipped
`tsd_search` and `tsd_report`, leaving no `3_TSD_search/pangenome.vcf`.

---

## Verify the installation

Without inputs the workflow prints its banner and stops with the input message, which shows
that Nextflow and the submodule are in order:

```bash
nextflow run cgroza/GraffiTE -r v1.1.0.0 -preview
```

```
V. 1.1.0.0 - v1.1.0.0
...
No input given. Pass one of --longreads, --bams, --assemblies, --pav, --svs (discovery), --vcf (a merged SV VCF), --RM_dir (RepeatMasker output of an earlier run) or --graffite_vcf (a pangenome.vcf from an earlier run).
```

<span class="src">`main.nf:151`</span>. Nextflow pulls the container only when a process runs,
so the quick install test below is the first check of the image.

### Quick install test

`test/quick_test.tar.gz` holds the inputs of a small run that uses the container: a 195 kb
synthetic genome in 8 contigs, 4 haploid assemblies and 2 short-read samples. The run goes through
discovery, annotation and genotyping, and `check.sh` compares the output with what GraffiTE 1.1
writes for this data set:

```bash
nextflow pull cgroza/GraffiTE -r v1.1.0.0
tar xzf ~/.nextflow/assets/cgroza/GraffiTE/test/quick_test.tar.gz && cd graffite_quick_test
nextflow run cgroza/GraffiTE -r v1.1.0.0 -params-file quick_test.yaml
bash check.sh
```

The last line it prints is `GraffiTE quick test: PASS`. Where Nextflow could use 4 CPUs and
12 GB of memory, the run took about 70 seconds once the image was downloaded. The test genotypes with
giraffe, because GraffiTE requests 32 GB for each task of PanGenie, the default back end.

---

## Updating

```bash
nextflow pull cgroza/GraffiTE -r v1.1.0.0
```

or add `-latest` to any `nextflow run` command. If Nextflow prints

```
cgroza/GraffiTE contains uncommitted changes -- cannot pull from repository
```

the cached copy was edited in place (usually `nextflow.config`). Remove it and pull again:

```bash
rm -rf ~/.nextflow/assets/cgroza/GraffiTE/
nextflow pull cgroza/GraffiTE -r v1.1.0.0
```

Each release names its own image tag in `nextflow.config`, so the first run of a new release
pulls the image built for it.

Nextflow does not check for a rebuilt image. Its Singularity cache is keyed on the
image **name**, so it reuses a `cgroza-graffite-1.1.0.0.img` left by an earlier run. A new
release names a new tag, and Nextflow pulls that one. If you point a run at `:latest` yourself,
the cached `cgroza-graffite-latest.img` stays whatever `:latest` was on the day you first pulled
it, and nothing in the log says the file is old. Delete it to force a fresh pull:

```bash
rm -f "${NXF_SINGULARITY_CACHEDIR:-$NXF_HOME/singularity}"/cgroza-graffite-latest.img
```

`NXF_HOME` is `~/.nextflow` unless you set it.

To see which image you are running, read the digest out of the cached file and compare
it with the registry:

```bash
apptainer inspect /path/to/cgroza-graffite-latest.img | grep base.digest
```

Do not compare checksums of the `.img` files themselves: Apptainer's
conversion from the Docker image is not byte-reproducible, so two conversions of one image have
different checksums.

---

## Execution profiles

`-profile standard` runs everything on the local machine and `-profile cluster` submits each
process to SLURM with the per-process CPU, memory and time parameters. `-profile cloud` submits
each process to AWS Batch (`process.executor = 'awsbatch'`) with Docker in place of Singularity,
because AWS Batch runs Docker images only <span class="src">`nextflow.config:7-29`</span>. It needs
your own AWS settings in a config passed with `-c`: at least a Batch queue (`process.queue`), a
region (`aws.region`) and a work directory on S3 (`-w s3://...`); see Nextflow's
[AWS Batch documentation](https://www.nextflow.io/docs/latest/aws.html). `pav_asm` cannot run
there, since its image is a Sylabs library image. We have not tested this profile.
[Resources and scaling](../guides/resources.md) covers how to size the requests.
