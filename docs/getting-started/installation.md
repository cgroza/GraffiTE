---
title: Installation
description: >-
  What GraffiTE needs on the machine, how to get the pipeline and its container
  image, and how to check the install before a real run.
---

# Installation

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `ee7da10`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](v1.0-vs-v1.1.md).

GraffiTE is a Nextflow pipeline. Every tool it calls lives in one container image, so the host
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

    Nextflow downloads and caches the pipeline the first time you run it, and initialises the
    `panmethyl` submodule for you:

    ```bash
    nextflow pull cgroza/GraffiTE -r v1.1dev
    ```

    The cached copy lives under `~/.nextflow/assets/cgroza/GraffiTE/`. Later runs use
    `nextflow run cgroza/GraffiTE -r v1.1dev ...`.

=== "Clone the repository"

    ```bash
    git clone --recurse-submodules https://github.com/cgroza/GraffiTE.git
    cd GraffiTE
    git switch v1.1dev
    ```

    `--recurse-submodules` matters. `main.nf` includes the `panmethyl` submodule unconditionally
    <span class="src">`main.nf:12`</span>, and a plain clone leaves that directory empty, so
    Nextflow stops at compile time with:

    ```
    Error main.nf:12:1: Invalid include source: '<dir>/panmethyl/module.nf'
    ```

    Run `git submodule update --init` in the clone and start again. `nextflow pull` clones
    fetch the submodule themselves.

---

## Get the container image

The config points every process at `docker://cgroza/graffite:latest`
<span class="src">`nextflow.config:10,15,21`</span>, and Nextflow pulls it on first use.
Apptainer converts the Docker image to a SIF itself, so nothing else has to be registered.

The image is `linux/amd64` only. On an arm64 machine it runs under emulation, slowly.

If the compute nodes have no internet access, pull the image once on a node that does and point
the run at the file:

```bash
apptainer pull graffite_latest.sif docker://cgroza/graffite:latest
nextflow run cgroza/GraffiTE -r v1.1dev -with-singularity /abs/path/graffite_latest.sif ...
```

`-with-singularity` overrides the image path in `nextflow.config`.

The PAV entry point (`--pav`) uses a second image, `library://becklab/pav/pav:latest`
<span class="src">`nextflow.config:339`</span>, pulled the same way. What the image contains,
and how the recipe in `GraffiTE.def` relates to it, is on [Container contents](../reference/container.md).

### Paths and `/tmp`

The config runs every container with `--contain --bind <dir>:/tmp`
<span class="src">`nextflow.config:175`</span>. Two consequences:

- The host filesystem is hidden from the container except for what Nextflow mounts: the work
  directory and the files it stages (`singularity.autoMounts = true`
  <span class="src">`nextflow.config:4`</span>). Give input files, including the paths inside
  samplesheets, as absolute paths.
- `/tmp` inside the container is the launch directory unless you say otherwise. Where that
  filesystem is small, slow, or not writable from inside the container, point `--container_tmp`
  at one that is:

  ```bash
  nextflow run cgroza/GraffiTE -r v1.1dev --container_tmp /scratch/you/tmp ...
  ```

  The directory must exist before the run starts. On SLURM sites that set a per-job scratch
  variable, `--container_tmp $SLURM_TMPDIR` works too.

!!! warning "Do not edit the cached config to do this"
    Before v1.1 the only route was to edit `singularity.runOptions` in
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
nextflow run cgroza/GraffiTE -r v1.1dev -preview
```

```
V. 1.1.0 - v1.1dev
...
No input given. Pass one of --longreads, --bams, --assemblies, --pav, --svs (discovery), --vcf (a merged SV VCF), --RM_dir (RepeatMasker output of an earlier run) or --graffite_vcf (a pangenome.vcf from an earlier run).
```

<span class="src">`main.nf:143`</span>. The container is only pulled when a process runs, so
the [Quickstart](quickstart.md) is the first check of the image.

---

## Updating

```bash
nextflow pull cgroza/GraffiTE -r v1.1dev
```

or add `-latest` to any `nextflow run` command. If Nextflow answers

```
cgroza/GraffiTE contains uncommitted changes -- cannot pull from repository
```

the cached copy was edited in place (usually `nextflow.config`). Remove it and pull again:

```bash
rm -rf ~/.nextflow/assets/cgroza/GraffiTE/
nextflow pull cgroza/GraffiTE -r v1.1dev
```

The image and the code are updated separately. After a code update that changes a tool, pull
the image again as well; the [Changelog](../changelog.md) says when that is needed.

---

## Execution profiles

`-profile standard` runs everything on the local machine and `-profile cluster` submits each
process to SLURM with the per-process CPU, memory and time parameters. `-profile cloud` sets
`process.executor = 'aws'` <span class="src">`nextflow.config:7-24`</span>, which Nextflow does
not recognise, so the run stops with `Unknown executor name: aws`. The AWS Batch executor is
called `awsbatch`, and needs your own AWS settings in a config passed with `-c`. All three
profiles name the same image. How to size the requests is on
[Resources and scaling](../guides/resources.md).
