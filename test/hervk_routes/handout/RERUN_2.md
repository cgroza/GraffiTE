# Re-run 2 — instructions

> **Superseded. Read `E2E_TEST.md` instead.** Re-run 2 completed and passed
> (see `HERVK_RUN_STATE.md`); this file is kept as the record of what it was
> testing. Two of its expectations have since changed — the three
> `REF_ARCH_CONFLICT` records it flags for analysis turned out to be tandem
> duplications and are now classified as such, so they raise no conflict at
> all. Do not run against this file's checklist.

Run 1 (commit `28109e7`, local) **passed all 41 assertions**. Its results are
still not usable, because analysing them turned up five bugs, three of which
change the output. This is the corrected re-run.

Read `README_HPC_CLAUDE.md` for the general setup; this file covers only what
is different.

## What changed since run 1

Pipeline is now at **`0681c33`** on `v1.1dev-hervk-v2`, pushed to GitHub. It
includes both of the patches you made locally in run 1 (the `tile_hits`
one-pass count and the `HERVK_REF_RM_OUT` resume hook), so `PATCHES.txt` is no
longer needed — please drop it.

| # | fix | changes output? |
|---|---|---|
| 1 | `ARCH_2LTR`/`ARCH_PERM` are polarity-aware. They were inverting `HERVK_ALLELE_REF`/`HERVK_ALLELE` on every `ARCH_*` **deletion**. | **yes** |
| 2 | `is_ltr_family` is an explicit HML-2 allowlist, not `startswith('LTR5')`. Your finding — `LTR57-int`, `LTR53-int`, `LTR54B` were being counted as LTR bp. | **yes** |
| 3 | Reference windows are re-cut when the element reaches a window edge. | **yes** |
| 4 | `|SVLEN|` cap of 25 kb on candidacy. Your finding — this is what let the 25.3 Mb deletion in. | no (drops junk) |
| 5 | `tile_hits` counts in one pass. Your fix, folded in. | no |

## Why the reference has to be masked again

Fix 3 is the reason, and it is worth being precise about, because the other
four do not need it.

At `flank=1500` an **insertion** footprint is a point, so its window is ~3 kb.
A complete HML-2 provirus is 9,472 bp. The window physically cannot contain the
element it is being asked to measure, so a proviral reference locus reads
`partial`. Run 1 has 22 of those.

That is not a parsing problem — the sequence outside the narrow window was
never masked, so there is nothing on disk to re-read. The fix cuts wider
windows (`flank=12000`) for exactly those candidates and masks *that* sequence.
This is why the rescue pass deliberately ignores `HERVK_REF_RM_OUT`: a `.out`
computed over the run-1 windows cannot describe the wider ones.

Concretely, this is what made chr6:78,894,316 look self-contradictory in run 1:

```
chr6-78894317-DEL-8465   ref_state=provirus   elem span 9425   window spans the deletion
chr6-78894876-INS-8465   ref_state=partial    elem span 2058   window is 3 kb, element is 9.5 kb
```

Same locus, same reference, opposite answers — purely an artefact of window
size.

**The other four fixes need no masking.** Fixes 1, 2 and 5 are code over tables
that already exist; fix 4 only removes candidates.

## Why this will be fast now

Run 1 sent 26.4 Mb to RepeatMasker and the masking phase alone ran 1h39m before
the 4h wall killed `ProcessRepeats`. 25.3 Mb of that — 95.8% — was one window,
around the 25.3 Mb deletion that fix 4 now excludes from candidacy.

Expected this time: roughly 1.1 Mb for the first pass, plus ~30 rescue windows
at ~24 kb each (~0.7 Mb). Under 2 Mb total, versus 26.4 Mb. Minutes, not hours.

So: **do not try to reuse anything.** Leave `HERVK_REF_RM_OUT` unset and let it
mask clean. The resume hook stays in the code for future runs, but using it
here would only reintroduce the stale narrow windows.

## Run it

```bash
cd /xdisk/cgoubert/cgoubert/GraffiTE1.1/CaG
./bootstrap.sh                  # pulls 0681c33, refreshes the handout
                                # (your INPUTS.env is preserved)
unset HERVK_REF_RM_OUT
./preflight.sh
./run_hervk_test.sh
./bundle_results.sh
```

`INPUTS.env` needs no changes — same `RM_DIR`, `REFERENCE`, `TE_LIBRARY`.

### -resume is handled for you

`run_hervk_test.sh` now stamps the pipeline commit into `$OUTDIR/.last_commit`
and **drops `-resume` automatically when the commit moves**. This matters more
than it looks: Nextflow's task hash covers the process script and the input
files but **not the contents of `bin/`**, which is staged onto `PATH`. The
`hervk_annotate` process block is unchanged between run 1 and run 2, so a plain
`-resume` would have hit the cache and replayed run 1's HERV-K output while
reporting success. You will see:

```
!! pipeline moved 28109e7 -> 0681c33
!! running WITHOUT -resume: ...
```

If you do not see that line and `$OUTDIR` already holds run-1 output, stop and
say so — it means the guard did not fire.

## What must be true this time

Everything from `README_HPC_CLAUDE.md` still applies, plus:

1. **Zero `REF_ARCH_CONFLICT` notes.** New `HERVK_NOTE` field in
   `hervk_calls.tsv`. It fires when the architecture and the masked reference
   imply different REF states. Run 1 (replayed with fix 1) had two; after the
   fix, zero. Any conflict remaining is a genuine disagreement and is the most
   interesting thing in the output — report it in full.

2. ~~**`chr6-78894876-INS-8465` should raise `REF_ARCH_CONFLICT`**~~ — no
   longer true. That record is a tandem duplication (a second proviral unit in
   an LTR of the existing provirus), is now classified `tandem_prov`, and has
   its genotypes withheld. Both chr6 records agree the reference is a provirus,
   so there is no conflict to raise. See `E2E_TEST.md`.

3. **`partial` count should drop sharply** from 22. Any `partial` that survives
   the rescue pass is a real partial element, not a windowing artefact.

4. **`null_prov` still 3**, `ARCH_PERM` k values unchanged (102, 174, 303, 410,
   574, 707, 735, 868), the same merge loci flagged.

5. **`other` should be lower** than run 1's 79 — fixes 2 and 3 both move records
   out of it. On the human subset the replay gives zero.

## Report back

The bundle, plus:

- the assertion log
- any `HERVK_NOTE` values that are not `.`
- the `ref_state` distribution (run 1: null=134, solo=45, partial=22, provirus=7)
- how many windows the rescue pass re-cut (stderr line from `hervk_ref_state`)
- wall time for `hervk_annotate` from `nextflow_trace.txt`

Do not tune parameters to make assertions pass. A disagreement is the result.
