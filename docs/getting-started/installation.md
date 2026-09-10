---
title: Installation
description: >-
  What GraffiTE needs on the machine, how to get the pipeline and its container
  image, and how to check the install before a real run.
---

# Installation

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
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

!!! warning "Nextflow 26 and the strict parser"
    Nextflow 26.04.6 rejects `main.nf` with its default parser: it stops at a line continuation
    written with a trailing dot (`Unexpected input: 'splitText'`) and at the `switch` statements.
    Set the legacy parser before running:

    ```bash
    export NXF_SYNTAX_PARSER=v1
    ```

    With that variable set, Nextflow 26.04.6 parses and previews the workflow. Older releases
    that still default to the legacy parser need nothing.

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
    <span class="src">`main.nf:42`</span>, and a plain clone leaves that directory empty. Since
    commit `959044b` the workflow stops with:

    ```
    panmethyl/ is empty: the submodule is not initialised. Run `git submodule update --init` in <dir>, or clone with --recurse-submodules.
    ```

    <span class="src">`main.nf:38-40`</span>. Run the command it names and start again.

---

## Get the container image

The config points every process at `library://cgroza/collection/graffite:latest`
<span class="src">`nextflow.config:11,16,22`</span>, and Nextflow pulls it on first use. The Sylabs
library must be registered with Apptainer once:

```bash
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud
```

If the compute nodes have no internet access, pull the image once on a node that does and point
the run at the file:

```bash
apptainer pull --arch amd64 graffite_latest.sif library://cgroza/collection/graffite:latest
nextflow run cgroza/GraffiTE -r v1.1dev -with-singularity /abs/path/graffite_latest.sif ...
```

`-with-singularity` overrides the image path in `nextflow.config`. The same image is on Docker
Hub and Apptainer can pull it from there:

```bash
apptainer pull graffite_latest.sif docker://cgroza/graffite
```

The PAV entry point (`--pav`) uses a second image, `library://becklab/pav/pav:latest`
<span class="src">`nextflow.config:322`</span>, pulled the same way. What the image contains,
and how the recipe in `GraffiTE.def` relates to it, is on [Container contents](../reference/container.md).

### Paths and `/tmp`

The config runs every container with `--contain --bind $(pwd):/tmp`
<span class="src">`nextflow.config:5`</span>. Two consequences:

- The host filesystem is hidden from the container except for what Nextflow mounts: the work
  directory and the files it stages (`singularity.autoMounts = true`
  <span class="src">`nextflow.config:4`</span>). Give input files, including the paths inside
  samplesheets, as absolute paths.
- `/tmp` inside the container is the launch directory. If that filesystem is small or slow,
  edit `singularity.runOptions` in your copy of `nextflow.config` (or in
  `~/.nextflow/assets/cgroza/GraffiTE/nextflow.config`) to bind a larger writable directory:

  ```groovy
  singularity.runOptions = '--contain -B /scratch/you/tmp:/tmp'
  ```

---

## Verify the installation

Without inputs the workflow prints its banner and stops with the input message, which shows
that Nextflow, the parser setting and the submodule are in order:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run cgroza/GraffiTE -r v1.1dev -preview
```

```
V. 1.1.0 - v1.1dev
...
No input given. Pass one of --longreads, --bams, --assemblies, --pav, --svs (discovery), --vcf (a merged SV VCF), --RM_dir (RepeatMasker output of an earlier run) or --graffite_vcf (a pangenome.vcf from an earlier run).
```

<span class="src">`main.nf:151`</span>. The container is only pulled when a process runs, so
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

`-profile standard` runs everything on the local machine, `-profile cluster` submits each
process to SLURM with the per-process CPU, memory and time parameters, and `-profile cloud`
targets AWS Batch <span class="src">`nextflow.config:8-25`</span>. All three use the same
image. How to size the requests is on [Resources and scaling](../guides/resources.md).
