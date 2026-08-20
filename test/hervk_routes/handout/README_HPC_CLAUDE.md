# HERV-K v2 discovery test — handout

You are running one test on a cluster and reporting what came back. Everything
needed is in this directory. Read this file top to bottom before starting.

## The one thing not to do

**Do not run graph genotyping.** There are no raw reads for this cohort. The
`--genotype false` flag in `run_hervk_test.sh` is deliberate. If a step appears
to want reads, something is wrong with the invocation — stop and report it
rather than looking for FASTQs.

## What is being tested

`bin/hervk_classify.py` was rewritten. The old version decided whether a HERV-K
locus was solo-vs-provirus or null-vs-provirus from the size of the SV alone.
That cannot separate an 8.5 kb LTR+internal block inserting into an existing
solo LTR from a complete provirus inserting into an empty site, and it made one
class (`null_prov`) unreachable.

The new version reads the raw RepeatMasker tables to recover the element's
internal architecture, and masks a window of the reference genome at each
candidate to see what the reference actually holds.

The expectations in `EXPECTED.tsv` were measured from the RepeatMasker tables
of the *existing* CaG run before the new code was written. So this is a real
test: the answers were derived independently, and the classifier is being asked
to reproduce them.

## Inputs

Everything is driven by **`INPUTS.env`** in this directory. Populate the three
paths at the top; the rest have working defaults.

| variable | what |
|---|---|
| `PAV_VCF` | the merged PAV call set for the CaG cohort |
| `REFERENCE` | CHM13v2 FASTA (a `.fai` beside it helps but is built if absent) |
| `TE_LIBRARY` | RepeatMasker library FASTA — **must contain `LTR5_Hs` and `HERVK-int`** |
| `GRAFFITE_SIF` | optional local `.sif`; leave empty to let Nextflow pull the image |
| `OUTDIR`, `PROFILE`, `CPUS` | defaults `hervk_v2_run`, `cluster`, `8` |
| `REVISION` | GraffiTE branch — `v1.1dev-hervk-v2` |

All paths must be **absolute**. `nextflow.config` runs the container with
`--contain`, so a relative path resolves to nothing inside it and the run fails
late, after the allocation is already spent. Pre-flight warns about this.

## Run it

The pipeline is not checked out locally — Nextflow fetches and caches
`cgroza/GraffiTE` at the pinned revision, and `bootstrap.sh` copies this
handout out of that cache.

```bash
cd /xdisk/cgoubert/cgoubert/GraffiTE1.1/CaG
module load nextflow            # whatever the site provides

./bootstrap.sh                  # nextflow pull + refresh this handout
$EDITOR INPUTS.env              # fill in the three paths
./preflight.sh                  # stop here if it fails
./run_hervk_test.sh             # pipeline, then assertions
./bundle_results.sh             # -> hervk_v2_results_<date>.tar.gz
```

`bootstrap.sh` never overwrites a populated `INPUTS.env` — a newer template
lands as `INPUTS.env.new` instead.

`run_hervk_test.sh` passes `-latest`, so a fix pushed upstream is picked up
rather than silently running a stale cache, and `-resume`, so a re-run after a
transient failure continues where it stopped.

Run it under `tmux`/`screen` or as a batch job: with `-profile cluster`
Nextflow submits its own slurm jobs and the driver process must stay alive.

Expected runtime: this annotates an existing call set and does not genotype, so
the long pole is RepeatMasker over the SV sequences — hours, not days. The
HERV-K step itself masks ~50 windows of ~12 kb and takes seconds.

## What must pass

`assert_hervk_test.py` exits non-zero on any failure. The two results that
decide whether this work is sound:

1. **`null_prov` is non-zero.** Specifically `chr19-21797327-INS-9478` and
   `chr8-7226885-INS-9468`, both previously called `solo_prov`. Their
   RepeatMasker architecture is `LTR[1-968] · INT[1-7536] · LTR[1-968]` — two
   whole LTRs, which only happens when nothing was consumed by the alignment,
   which means the reference was empty. The first of the two is independently
   labelled `pro_pre` (provirus vs pre-integration site) by Wildschutte 2016.

   If this comes back 0, the architecture layer is not reaching the classifier.
   Check that `3_TSD_search/hervk_arch.tsv` has rows.

2. **Three loci flagged `MERGE_CANDIDATE`:** `HERVK_chr11_101704640`,
   `HERVK_chr12_55299985`, `HERVK_chr6_78894316`. The chr11 one matters most —
   two records 574 bp apart that describe one insertion, which the existing
   snarl-based grouping misses because they sit in different snarls.

Also checked: eight `ARCH_PERM` records must reproduce their exact permutation
points (`k` = 102, 174, 303, 410, 574, 707, 735, 868), 22 lone LTRs must come
back `null_solo`, and the degenerate records must reach `REF_ANNOT`.

## Reading a failure

| symptom | what it means |
|---|---|
| `null_prov count is 0` | architecture not reaching the classifier — is `hervk_arch.tsv` populated? Are the `repeatmasker_dir` inputs staged? |
| many `evidence=UNRESOLVED` | same cause, or the RepeatMasker tables were not passed in |
| `ref_state=unknown` on `REF_ANNOT` records | reference masking produced nothing — check `TE_LIBRARY` really contains `LTR5_Hs`/`HERVK-int`, and that `samtools faidx` worked |
| `k` differs from expected | **report it, do not adjust anything.** `k` is an alignment property and is expected to vary between callers. If it moved under the same PAV input, that is a real finding |
| `HERVK_chr11_101704640` not flagged | locus grouping is not seeing both records — check `hervk_calls.tsv` contains both `chr11-101704641-INS-8504` and `chr11-101705464-INS-8498` |
| `POLARITY_CONFLICT` on chr6 | **expected.** The two chr6 records genuinely disagree about the reference state; that is what the flag is for |

## What to report back

Send the bundle plus a short summary:

- the assertion log verbatim (pass or fail)
- the `null_prov` records found
- the `MERGE_CANDIDATE` loci and their flags
- the class counts from `hervk_polymorphism_summary.md`
- **the reference state called at `chr6:78,894,316`** — this one is an open
  question, not a test. Two records there imply opposite reference states, and
  the masked reference is what settles it. Report `ref_state` for
  `chr6-78894317-DEL-8465` from `hervk_refstate.tsv` whatever it says.
- anything in the `NOTE:` lines of the assertion output

Do not tune parameters to make assertions pass. If something disagrees, that is
the result.
