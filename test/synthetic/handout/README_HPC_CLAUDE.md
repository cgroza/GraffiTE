# Synthetic end-to-end test set — handout

You are building and running the first end-to-end validation of GraffiTE v1.1 on
a cluster, and reporting what came back. Everything needed is in this directory.
Read this file top to bottom before starting.

Unlike the earlier handouts (`test/pangenie_93/`, `test/hervk_routes/`), there is
no reporter's data to point at. **This handout generates its own inputs.** You
fill in where things live, not what they are.

## Why this exists

v1.1 is 369 commits ahead of `main` and has never been run end to end. Every
feature has a unit test over fixtures; nothing has driven discovery →
annotation → TSD → graph → genotyping in one pass. The release PR
(cgroza/GraffiTE#105) was withdrawn for that reason and is waiting on this.

The second goal is documentation. The reference pages were written from the code
— `docs/getting-started/quickstart.md:75` says so in as many words — so every
claim about outputs is a claim about what the code appears to do. The runs here
are what lets those pages be corrected against evidence.

The full plan is `docs/design-notes/synthetic-test-set.md` in the repo. This
handout is Tier 1 and Tier 2 of it.

## The one thing not to do

**Do not tune anything to make an assertion pass.** If a number disagrees, that
is the result. Several expectations here are deliberately *not* frozen yet,
because they depend on tool heuristics nobody controls — see
[Calibration](#calibration-comes-first).

## What gets built

`preflight.sh` runs `build_synthetic.py`, which writes `$WORKDIR/build/`:

| | |
|---|---|
| `ref/synth.fa` | 8 contigs, ~200 kb. `chr1` is the annotation zone, `chr2` the HML-2 zone, `chr3_quiet` has reference TEs and no variant, `chr4_narrow` has one symbolic `<INV>` and no indel, and `X`/`chrX`/`chrX_alt`/`chrY` carry the same element for the ploidy probe. |
| `lib/synth_TE.fasta` | 7 entries. Consensus from the `human_DFAM3.6.fasta` inside `test/human_test_set.tar.gz`, except L1HS: that library ships `L1HS_5end` and `L1HS_3end` separately and no full-length copy, so this one joins them. |
| `lib/synth_TE_bare_HERVK.fasta` | Identical but for the internal region, named `HERVK` rather than `HERVK-int`. See [the dual-library cell](#the-dual-library-cell). |
| `hap/h{1..4}.fa` | Four haploid assemblies. Every planted site is carried by two of the four, so nothing is fixed. |
| `reads/`, `bam/` | Short and long reads off the haplotypes, and two coordinate-sorted BAMs with `@RG SM` set. Error-free on purpose: a genotyping failure should be the graph's, not the read set's. |
| `vcf/` | Hand-written fixtures for `--svs` and `--vcf`, carrying the real inserted sequence so RepeatMasker has something to annotate. |
| `truth.tsv` | One row per planted site: contig, position, what it is for, length, TSD. |
| `MANIFEST.sha256` | Content hashes. Samplesheets and BAMs embed the absolute output path and are listed separately as path-dependent. |

The generator is deterministic: same `SEED` and same Dfam source, same bytes.
Verified by building twice into different directories and diffing the manifest.

## Inputs

Everything is driven by **`INPUTS.env`**. The fields that matter:

| variable | what |
|---|---|
| `WORKDIR` | absolute, on `/xdisk`. `work/` will outgrow a home quota. |
| `CONTAINER_TMP` | absolute, on scratch. See the warning below. |
| `GRAFFITE_SIF` | optional. Worth setting: the image tag is mutable. |
| `PROFILE` | `standard`. Not `cluster`. See below. |
| `CPUS`, `MEM_GB` | must match `--cpus-per-task` and `--mem` in `submit_driver.sh`. |
| `RUNS` | which matrix cells to run by default. |

All paths absolute: `nextflow.config` runs the container with `--contain`, so a
relative path resolves to nothing inside it and the run fails late, after the
allocation is spent.

!!! Two things that have bitten this project before

**`CONTAINER_TMP`.** If `/tmp` inside the container is not writable or fills,
the run does not fail. It finishes having skipped `tsd_search` and `tsd_report`,
with no `3_TSD_search/pangenome.vcf`, and reports success. `--container_tmp` is
new in v1.1; before it the only route was editing the cached `nextflow.config`,
which then makes `nextflow pull` and `-latest` fail with `contains uncommitted
changes`.

**`PROFILE=standard`, inside `sbatch`.** Not `cluster`. On Puma the cluster
profile submits one slurm job per task and has cost roughly 85 minutes of queue
latency *each*; job 23613713 hit the 12 h wall having finished a fraction of the
work. The pattern that works here is one allocation, the local executor, and
`local.config` sized to the allocation. The login node has no `apptainer`, so
`preflight.sh` will fail there — run it inside the allocation, or accept that it
reports the missing runtime and fix that one item later.

## Run it

```bash
mkdir -p /xdisk/cgoubert/cgoubert/GraffiTE1.1/synthetic && cd $_
module load nextflow

nextflow pull cgroza/GraffiTE -r test/synthetic-end-to-end
cp ~/.nextflow/assets/cgroza/GraffiTE/test/synthetic/handout/* .

./bootstrap.sh                  # pull + refresh this handout
$EDITOR INPUTS.env              # WORKDIR, CONTAINER_TMP, sizes
./preflight.sh                  # checks, then builds the inputs
sbatch submit_driver.sh         # the whole matrix
# or, inside an allocation:
./run_matrix.sh spine
```

`./run_matrix.sh -l` lists the cells. Each is independent except that everything
after `spine` reuses the spine's published outputs, so run `spine` first and let
it finish. A failing cell does not stop the others; per-cell status goes to
`MATRIX.log` and each cell's own log to `$WORKDIR/runs/<cell>/run.log`.

## Calibration comes first

**Before freezing any expectation, run `spine` and read what it produced.**

The assertions that matter rest on tool behaviour nobody here controls: whether
RepeatMasker reports a planted copy at a given divergence, how ProcessRepeats
assigns link IDs, whether truvari collapse merges two HERV-K records ~968 bp
apart, whether PanGenie's k-mer model fits a genome three orders of magnitude
smaller than it was written for.

So the first pass is a **measurement**, not a test. From the spine's output,
record into `OBSERVED.tsv`:

- every planted site in `truth.tsv`: did it reach `pangenome.vcf`, with what
  `n_hits`, `repeat_ids`, `matching_classes`, `total_repeat_span`, `TSD`,
  `polyA`, and what SW score in the RepeatMasker `.out`
- the `L1_5PINV` locus (`A09_l1_twin_primed_Cplus`): the link IDs its fragments
  got, and whether `L1_5PINV` is set. **Record `None` if it is not.** The
  library has one full-length L1HS entry precisely so this can work, but
  ProcessRepeats decides, not us.
- the SVA VNTR-only locus (`A12_sva_vntr_only`): the reported
  `target_start`/`target_end`, and whether they land strictly inside SVA_E's
  429..863 window
- the span ladder (`A14`..`A21`): which rungs survived `--repeat_span_cutoff`

Then check these, and treat a failure as a **fixture** problem with its own
message, not a pipeline regression:

- every designed element produced at least one RepeatMasker hit
- no planted copy sits within 15 % of the lowest observed SW score
- the span ladder has exactly one cut point, and the rungs either side of it
  bracket 0.80
- `hervk_loci.tsv` has at least one row with `n_records > 1` — the whole
  consolidation layer depends on truvari *not* merging the designed multi-record
  locus, and nothing else guards that

If calibration fails, say so and stop. Regenerating the genome with different
divergences is the fix, and it has to happen before 31 loci are committed to.

## What the matrix is for

| cell | what it is the only cover for |
|---|---|
| `spine` | discovery from assemblies, the truvari collapse path (≥2 VCFs, which is the only path that runs `bcftools norm -f`), RepeatMasker, TSD, polyA, `--human`, the whole HERV-K stack, giraffe genotyping |
| `pangenie` | the `-N` left-alignment guard added in `a106944`, and non-zero `not_in_graph`/duplicate counts in `genotyping_record_audit.tsv` |
| `graphaligner` | `vg construct` + GraphAligner |
| `precomputed` | `--graph` reuse, building nothing |
| `longreads` | two minimap2 presets, `sniffles` per sample then the joint call |
| `bams` | the no-`type`-column entry, and `@RG SM` deciding sample names |
| `vcf` | the `from_vcf=true` branch: no `tabix`, INFO kept, `SVLEN` **not** recomputed |
| `svs` | the other single-input branch, which *does* `tabix` |
| `duallib` | see below |
| `guards` | three launch-time errors nothing else tests |

`chr4_narrow` is worth knowing about: `bcftools index -s` lists only contigs that
carry records, so a contig with one symbolic `<INV>` and no indel *does* become a
chunk, `bcftools view --types indels` returns nothing for it, and the no-indel
exit in `bin/repmask_vcf.sh` fires. That is the only end-to-end route to the
branch whose absence crashed `tsd_prep` in issue #93.

### The dual-library cell

The Dfam library this repo ships names the HERV-K internal region **`HERVK`**.
The `--human` pair rule at `module/main.nf:548` requires `repeat_ids~"^HERVK-int"`.
So a user following the quickstart with the shipped library gets no pair-rule
records at all, silently.

`duallib` runs the same genome against `synth_TE_bare_HERVK.fasta` and the
expectation is that the carve-out records **vanish with no error and no
warning**. If they do, that is a real finding about what users get, not a test
failure. Report it either way.

## Reading a failure

| symptom | what it means |
|---|---|
| run finishes, no `3_TSD_search/pangenome.vcf` | `/tmp` in the container. Set `CONTAINER_TMP` and rerun. Check `runs/<cell>/run.log` for whether `tsd_search` ran at all. |
| `cp: cannot stat 'repeatmasker_dir/repeatmasker_dir/*'` | the issue #93 crash is back; `6864e30` should have removed that line |
| Nextflow exits in ~3 s with a config error | usually `-latest` against a dirty cached checkout. `preflight.sh` warns about this. |
| scheduler kills the driver | `CPUS`/`MEM_GB` in `INPUTS.env` disagree with `submit_driver.sh`, so the local executor oversubscribed the allocation |
| PanGenie runs but genotypes nothing | the most likely real failure here. It is written for gigabase genomes and GraffiTE exposes no `k`. Report the genotype counts; do not tune. |
| `vg autoindex` refuses the VCF | contig order. The generator writes them in ASCII-lexicographic order for this reason; if you edited the contig list, that is why. |
| `hervk_reconcile` refuses the back end | expected on anything but giraffe. It is a guard, not a bug. |

## What to report back

`./bundle_results.sh` packs it. Send the archive plus a short summary:

- `MATRIX.log`: which cells ran, which failed
- `OBSERVED.tsv` and the calibration verdict, **especially anything recorded as
  `None`** — those are the checks that were designed in and did not fire
- `genotyping_record_audit.tsv` from the `pangenie` cell: the `lost_at` counts.
  Non-zero `not_in_graph` is expected; zero means the doctored records did not
  make it in
- the `duallib` result, in the terms above
- `VERSIONS.txt` — `docs/reference/container.md` currently *guesses* at the tool
  versions and this is what replaces the guesses
- `published.txt` per cell; these become the outputs-documentation check
- wall time and peak memory per stage from `nextflow_trace.txt`

Anything that disagrees with `docs/reference/outputs.md` or
`docs/reference/vcf-fields.md` is the point of the exercise, not a problem with
the run. Quote it.
