---
title: Synthetic end-to-end test set
description: >-
  A plan for a synthetic genome and the runs over it that validate every entry point,
  every graph method and --human before v1.1 is merged into main.
---

# Synthetic end-to-end test set

!!! warning "Proposed. Nothing here is built."
    This describes work to be done, not behaviour GraffiTE has. Status of each piece is
    in [Sequencing](#sequencing).

## Context

v1.1 is 369 commits ahead of `main` and has never been run end to end as a release
candidate. Every feature has a regression test over fixtures, and the two HPC handouts
under `test/` run the real stack on real genomes, but nothing in between exists: no
single input set that drives discovery, annotation, TSD, the graph and genotyping in one
pass with outputs we can diff. The release PR was withdrawn for that reason.

The other half of the problem is that the documentation was written from the code. Every
factual claim about outputs is a claim about what the code appears to do, not about what a
run produced. `docs/getting-started/quickstart.md:75` says so in as many words. A test set
that produces real outputs is what lets those pages be corrected against evidence.

### Constraints, measured

| | |
|---|---|
| Container | `docker://cgroza/graffite:latest`, `linux/amd64` only, single layer, 2.45 GB compressed, config created 2025-12-08. A mutable tag: the image can change with no commit here. |
| This machine | `darwin/arm64`, no Apptainer, no Singularity, Nextflow not on `PATH`. The full stack runs under emulation or on Linux x86_64. |
| Read simulation | Neither `GraffiTE.def` nor `Dockerfile` installs a simulator. Reads must come from a generator we write. |
| RepeatMasker | `bin/repmask_vcf.sh:75` passes `-lib` only, no Dfam and no `-species`, so annotation cost scales with the test library, not with a reference database. |

The arm64 point sets the shape of everything below: **only Tier 0 is a pre-push gate on
this machine.** Any claim that the full matrix runs "in minutes on a laptop" is a claim
about a machine nobody on this project has.

---

## Three tiers, three lifetimes

Merging these into one suite is what makes runtime estimates collapse.

### Tier 0: the pre-push gate

Container-free, no Nextflow, seconds. The twelve existing `test/*.sh` scripts plus four
additions. This is the only tier that runs natively here.

- The runner **counts executed checks and fails on zero**. Every current script exits 0
  on `[skip]`, so a machine missing `bcftools` reports success today.
- Fix `test/human_filter/run_test.sh:49`. Its pair clause is `n_hits==2`; the code at
  `module/main.nf:548` is `n_hits<=${params.hervk_pair_max_hits}`, default 3. The test
  mirrors a superseded expression and cannot catch a change to the threshold. Read the
  value from `nextflow.config` as the rest of that script already does, and add a 3-hit
  record to the fixture.
- Add a direct test of `hervk_reconcile`'s genotyper guard (`SUPPORTED_GENOTYPERS`).
- Add an N-in-ALT record to `test/pangenie_index/`'s fixture.
- Add a static check that every boolean parameter test in `main.nf` and `module/main.nf`
  routes through `isOn()`. That is the invariant commit `cfaff1e` established, and
  nothing enforces it.

### Tier 1: the spine

One expensive run pays discovery, RepeatMasker, TSD, polyA, the `--human` filter and the
HERV-K stack once. Three near-free genotyping re-entries follow it through
`--graffite_vcf` and `--graph`. This works because the whole discovery and annotation
block sits inside `if(!params.graffite_vcf)` at `main.nf:121`.

### Tier 2: the branch sweep

The flags no spine reaches, each asserting one thing: `--epigenomes` via a hand-written
`--lifted` CSV (which needs no methylation simulator), `--hervk_ref_annotation` re-run
against the spine's own reference `.out`, `--aligner winnowmap`, `--break_scaffolds` on an
N-containing haplotype, `--tsd_win 40`, and the `isOn()` negatives.

---

## The genome

Roughly 272 kb, plain uncompressed FASTA, under `test/synthetic/`. Length is nearly free, because
RepeatMasker never sees the genome, only `indels.fa` and the HERV-K reference windows. So
flanks are generous and the number of RepeatMasker invocations is not.

Contig names in ASCII-lexicographic order, which `vg autoindex` requires.

| Contig | Size | Purpose |
|---|---|---|
| `chr1` | 60 kb | The annotation zone. ~22 planted sites on a 2.5 kb pitch: TSD present and absent, polyA above and below `add_polyA.py`'s `MIN_LEN`, a reverse-complement polyT case, an `AluSx` that passes trusted and fails `--human`, the C+ twin-primed L1 and its `+C` negative, an SVA VNTR-only copy, and an eight-rung `total_repeat_span` ladder straddling 0.80. |
| `chr2` | 60 kb | The HML-2 zone. Six loci giving one reference state each: empty, solo LTR, full provirus, two-unit tandem array, truncated internal region, and a `HERVK9-int` decoy that must classify as non-HML-2. A seventh sits near the contig end with only 4 kb of flank, to exercise `hervk_ref_state.py`'s rescue pass and `truncated_by_window` on purpose. |
| `chr3_quiet` | 20 kb | Reference TE copies, no variant. Pins the per-chunk directory count. It does **not** reach the empty-chunk exit. |
| `chr4_narrow` | 20 kb | One symbolic `<INV>` record and no indel. `bcftools index -s` lists only contigs carrying records, so this contig *does* yield a chunk, `bcftools view --types indels` returns nothing, and `bin/repmask_vcf.sh`'s no-indel exit fires. This is the only end-to-end route to that branch. The fixture must declare every INFO tag it uses or bcftools exits 255. |
| `X`, `chrX`, `chrX_alt`, `chrY` | 4 kb each | The same planted `AluY` at the same offset, at 4/7/9 % divergence so alignments stay separable against `min_support = '2,4'`. `vg call -R chrX:1,chrY:1` is a full-match ploidy regex, not a name list: `chrX` and `chrY` must come back haploid, `X` and `chrX_alt` diploid. |

### The library is a hybrid, and its constants come from the code

`bin/hervk_arch.py` hard-codes `INT_CONSENSUS_LEN` and `LTR_CONSENSUS_LEN`. The shipped
`test/human_test_set.tar.gz` carries `HERVK#LTR/ERVK` at exactly that internal length and
`LTR5_Hs` at exactly that LTR length. Take HERV-K, SVA and Alu consensus from there.

L1 is different: that library has `L1HS_5end` and `L1HS_3end` as separate entries and no
full-length L1HS, so two library entries are likely to receive separate RepeatMasker
link IDs, which would make `L1_5PINV` unreachable. Calibration decides it. Hand-write a full-length L1HS for that locus
alone.

**The generator parses the two constants out of `bin/hervk_arch.py` and dies on
disagreement** rather than restating them in a comment.

### One arithmetic correction

An earlier draft claimed the SVA_D VNTR window admits no sub-copy that is both inside it
and above the 250 bp floor. That is wrong. `bin/annotate_vcf.R:133-134` gives SVA_D
`VNTR_start = 432`, `VNTR_end = 689`, and the test at `:144` is strict on both sides, so
the admissible band is 433..688, which is 256 bp. It clears the floor by six. Use **SVA_E**
anyway, whose 429..863 band is 435 bp wide, because RepeatMasker and not the generator
decides the reported hit edges and a six-base margin is not a margin.

---

## Calibrate before freezing anything

The assertions that matter rest on tool heuristics nobody here controls: whether
RepeatMasker reports a copy at a given divergence, how ProcessRepeats assigns link IDs,
whether truvari collapse merges two HERV-K records 968 bp apart, whether PanGenie's k-mer
model fits a genome three orders of magnitude smaller than it was written for.

So a calibration run **gates** the design rather than blessing it afterwards. Before any
expectation is committed, assert:

- every designed element produces at least one RepeatMasker hit
- no planted copy sits within 15 % of the lowest observed SW score
- the SVA_E VNTR-only copy's reported `target_start`/`target_end` land strictly inside the
  window
- the C+ L1 locus's fragments carry one link ID. **Recorded, not asserted**: write `None`
  down if they do not

A blessed wrong expectation is indistinguishable from a right one. If calibration fails,
the genome is regenerated before anything downstream is written.

Two exit codes: **1** for a broken structural invariant, **2** for "tool heuristics moved,
review and re-freeze". A re-freeze must be committed as a reviewable diff, so it appears
in `git log` rather than in someone's shell history. With a mutable `:latest` tag, exit 2
can fire on a rebuild nobody here initiated.

### Assert the fixture before asserting the pipeline

Some assertions test our own construction, not GraffiTE, and must fail with their own
message: that `hervk_loci.tsv` has a row with `n_records>1` (the whole consolidation layer
hangs on truvari not merging the designed multi-record locus), that the span ladder has
exactly one cut point, and that the two rungs nearest it straddle 0.80. Otherwise a
fixture that drifted reads as a classifier regression.

---

## The matrix

Discovery and annotation are paid once; everything after is cheap.

| Run | Covers |
|---|---|
| Spine | `--assemblies` ×4 + `--svs` + `--human` + giraffe + genotyping. Publishes the graph, the alignments, the `RM_dir` and the reference `.out` for reuse. |
| pangenie | `--graffite_vcf` on a doctored `pangenome.vcf` carrying N-in-ALT, overlapping and duplicate records, so the audit's `not_in_graph` and duplicate counts are non-zero. Pangenie-only: appending records leaves the file unsorted and `vg autoindex` would refuse it. |
| graphaligner | `--graffite_vcf`, long reads. |
| precomputed | Reuses the spine's graph and alignments. Both shapes: `--graph_alignments` and `--vcfs`. |
| longreads + bams | Two entry points, `sniffles_sample_call` then the joint call. A failure in either surfaces identically; split to diagnose. |
| vcf | The `from_vcf=true` branch: `cp`/`gunzip`, no `tabix`, INFO preserved, `SVLEN` **not** recomputed. |
| svs alone | Single-VCF, `from_vcf=false`: the branch that does `tabix` its input. |
| dual library | The same genome against a library spelling the internal region `HERVK` rather than `HERVK-int`, asserting the pair-rule records vanish from `pangenome.human.vcf` with no error and no warning. The repo's own shipped library uses the bypassing spelling, so this is what a user following the quickstart gets. |
| launch guards | `main.nf:46-49`, `:51-56`, `:61-63`. Deterministic, container-free, currently untested. |

Two vacuity guards on every run that reaches them, as hard assertions with their own
messages: `pangenome.vcf` record count above zero, and `pangenome.human.vcf` record count
above zero. A matrix of green existence checks over empty files is the failure mode this
whole exercise exists to avoid.

One judgement call worth recording: the giraffe "Falling back on single-end mapping"
warning is a property of genome size, not a GraffiTE defect. It is a **note** that triggers
deepening the read set, never a hard failure.

---

## Documentation

**Land the script half first, on its own schedule.** It needs no container, no Nextflow
and no synthetic genome, and it is the highest value per minute in this plan.

The largest class of documentation error is invisible to any run. Citations of the form
`module/main.nf:831` go stale whenever code moves, and `stamp_pages.py --set` rewrites the
"verified against commit X" stamp blind, so a page can carry a fresh stamp and wrong line
numbers at once. That already happened: all 23 process citations were correct at `67900f6`
and all 25 were wrong after PR #103 merged, off by 12 to 24 lines, with nothing to catch
it. `docs/scripts/remap_citations.py` fixed 343 of them; the repair is commit `9adea98`.

So:

1. Build the semantic checker as the **validating half of `remap_citations.py`**, sharing
   its `CITE` regex and skip list, rather than as a fourth script with a fourth regex. It
   asserts a cited line still holds what the row claims: that `module/main.nf:N` opens the
   named process, that `nextflow.config:N` assigns the named parameter. Both checks are
   already written ad hoc in `9adea98`'s verification and should be made permanent.
2. Guard `stamp_pages.py --set` against stamping a page whose citations are bad. Ten lines,
   and it stops a green stamp from being evidence of nothing.
3. Add a static `publishDir`-to-page check reading `module/main.nf` and
   `panmethyl/module/main.nf`. `docs/reference/outputs.md:3` calls itself "every file
   GraffiTE publishes" and lists no methylation output, while
   `docs/reference/processes.md` documents four. Note also that `check_params.py`'s
   `FIELD_SOURCES` omits `panmethyl/`, so `PMN`/`PML`/`PMD` were invisible to it while it
   reported ok.
4. Run all three advisory, clear the backlog, then make them blocking in
   `.github/workflows/docs.yml`.

**Then** the run half. Every edit names the run that justifies it:

- Delete `quickstart.md:75` and give it record counts from the spine.
- Correct the output-tree entries the runs refute, and add the four `4_Genotyping` files
  that v1.1 added.
- Add `precomputed` to the "giraffe is the only validated back end" warning if the
  precomputed run supports it.
- Replace `container.md`'s admitted version guesses with a `VERSIONS.txt` captured from the
  image.
- Stamp the pages **last**, after the citation checker passes.

The "documented but never produced" direction of the outputs check must be a warning keyed
to the flags the matrix ran, printing the unreached flag set. Otherwise
uncovered reads as absent.

---

## Sequencing

Each step is independently checkable; nothing later depends on an unverified earlier one.

1. Correct the stale premises in writing. Done in this document.
2. Land the documentation script half: the citation validator, the `publishDir` check, the
   `isOn()` check. Advisory, then blocking.
3. Guard `stamp_pages.py --set` on citation validity.
4. Fix and extend Tier 0. No new inputs.
5. Write the generator: hybrid library, constants parsed from `bin/hervk_arch.py`, the
   reference, four haplotypes, the read sets, the BAMs, the hand-written VCF fixtures.
   A second run must reproduce the manifest byte for byte, and no ID may exceed 50
   characters.
6. Run the calibration gate. If it fails, regenerate and repeat before writing anything.
7. Build and run the spine. Lift the per-sample `vg call` VCFs out of `work/` as an
   explicit numbered step, before any cleanup.
8. Assert the fixture, then the pipeline.
9. Run the three genotyping re-entries.
10. Run the remaining discovery entry points and the launch guards.
11. Run the dual-library falsification.
12. Run the Tier 2 branch sweep.
13. Feed every manifest to the outputs check and close the methylation gap.
14. Rewrite the documentation from the runs. Stamp last.
15. Reopen the `v1.1dev` → `main` PR, citing the matrix.

---

## Known risks

Things the plan depends on that nobody has verified. Each is a reason a green run might
mean less than it appears to.

- **PanGenie's k-mer model on a small genome.** Written for gigabases; GraffiTE exposes no
  `k` or coverage parameter; the container builds it from git master with no tag. If it
  genotypes nothing, every audit number the matrix quotes is meaningless. Check the
  genotype counts before trusting anything downstream.
- **Truvari collapse and the multi-record HERV-K locus.** Run with `--chain -P 0.5 -p 0.5
  -S -1 -k common`. Nobody checked its `refdist` against the installed binary. Guarded by
  the fixture assertion above, not by design.
- **ProcessRepeats link IDs** on a library this small. Both `L1_5PINV` and the whole HERV-K
  architecture layer depend on them.
- **RepeatMasker's default reporting threshold.** No `-cutoff`, `-div` or `-species` is
  passed, so whether a planted divergence is reported is a tool default nothing here pins.
- **`vg` under emulation** on arm64. If the release binary does not run, Tier 1 has no home
  on this machine at all.
- **The mutable `:latest` tag.** Only ULTRA, vg, pypy and Dfam are version-pinned inside
  the image.
- **`ARCH_PERM` may be unreachable.** It arises only when minimap2 and svim-asm place a
  provirus insertion inside a reference solo LTR, which no design can force.
- **Python's PRNG stability** across versions and platforms, which a committed manifest
  checksum is the only guard against.
