# Cluster test: the `--human` pME filter (branch `v1.1dev-human-filter`)

You are running a validation test of a change to GraffiTE's `--human` mode on real data.
This document assumes no prior context. Read it end to end before running anything.

---

## 1. What changed, and the one thing that can silently break it

`--human` used to take `pangenome.trusted.vcf` (a species-agnostic heuristic subset) and
filter it by TE **class**. It now builds its output **directly from `pangenome.vcf`** with
a filter written for human polymorphic mobile elements (pMEs), and keeps only young
subfamilies:

| class | kept | dropped |
|---|---|---|
| `SINE/Alu` | `AluY*` | `AluS*`, `AluJ*` |
| `LINE/L1` | `L1HS` | `L1PA*`, `L1M*`, `L1P*` |
| `Retroposon/SVA` | `SVA_D/E/F` | `SVA_A/B/C` |
| `Simple_repeat` | SVA-VNTR records (`SVA_[DEF](VNTR_only)`) | unrelated tandem repeats |
| `LTR/ERVK` | strict HML-2: `HERVK-int`, `LTR5_Hs`, `LTR5A`, `LTR5B` | `HERVK9/11/22/C4-int`, `MER11*`, `LTR13`, bare `LTR5` |

Two consequences to expect in the output:

- The file is named **`pangenome.human.vcf`** (was `pangenome.trusted.human.vcf`), and
  **`pangenome.trusted.vcf` is not produced at all** in `--human` runs.
- Records with `n_hits==2` annotated `HERVK-int` + SVA are now admitted. These are HML-2
  proviral SVs where RepeatMasker assigns part of the LTR to SVA (SVA/`LTR5_Hs`
  homology); the old `n_hits==1` rule discarded every one of them.

**The risk that makes §3 mandatory.** The filter is a single `bcftools view -i`
expression, and it depends on bcftools-specific regex behavior:

- `~` is evaluated **element-wise** on `Number=.` fields, with `^` anchoring **per
  element** — `repeat_ids~"^LTR5_Hs"` must match the record `SVA_A,LTR5_Hs`.
- bcftools regexes have **no alternation** (`|` and `\|` do not work inside a regex), which
  is why each whitelist is a comma-separated list expanded into OR'd clauses.

This was verified on bcftools 1.22. The GraffiTE container builds bcftools from
**unpinned git master** (`GraffiTE.def:195`), so the image may carry a different build.
If its semantics differ, the filter silently keeps the wrong records — it does not error.
§3 settles this in one command.

---

## 2. Inputs to locate, then confirm with the user

You need three paths. **Do not guess — resolve them, then echo all three back to the user
and wait for confirmation before launching.**

1. **`RM_dir`** — a previous GraffiTE run's `2_Repeat_Filtering/` directory. This is what
   lets the run skip SV discovery and RepeatMasker entirely.

   ```bash
   find <likely-parent-dirs> -maxdepth 4 -type d -name 2_Repeat_Filtering 2>/dev/null
   ```

   Verify the layout — it must contain numbered subdirectories, each holding both files:

   ```bash
   ls <RM_dir>                                        # -> 1  2  3  ...
   ls <RM_dir>/1/genotypes_repmasked_filtered.vcf     # must exist
   ls -d <RM_dir>/1/repeatmasker_dir                  # must exist
   ```

   If either is missing the pipeline fails at channel creation with a `checkIfExists`
   error.

2. **`--reference`** — the **same** reference FASTA used by the original run. `fix_vcf.py`
   pulls REF alleles from it; a different reference produces wrong REF fields.

3. **Previous human VCF (optional, but do look for it)** — the earlier run's
   `3_TSD_search/pangenome.trusted.human.vcf`. It is the before-picture for the
   comparison in §5g. If it does not exist, say so and skip that check.

---

## 3. Pre-flight: verify filter semantics inside the container — BLOCKING

```bash
git clone -b v1.1dev-human-filter https://github.com/cgroza/GraffiTE.git
# reuse the cached image if $NXF_SINGULARITY_CACHEDIR already has it, otherwise:
singularity pull graffite.sif library://cgroza/collection/graffite:latest
singularity exec graffite.sif bash GraffiTE/test/human_filter/run_test.sh
```

Expected — 7 `ok` lines and `all human filter tests passed`:

```
ok   defaults
ok   human_l1_ids adds L1PA2
ok   human_min_svlen=100 admits the 140 bp AluY
ok   hervk_sva_pair=false drops the HERVK+SVA pairs
ok   human_ignore_filter=true admits the non-PASS record
ok   human_hervk_ids=^LTR5 admits the ancestral LTR5
ok   empty human_sva_ids keeps the whole SVA/Simple_repeat classes
```

**Any `FAIL` here means the container's bcftools does not behave as the filter assumes.
Stop. Do not launch the pipeline.** Report which case failed, its expected-vs-got line,
and `singularity exec graffite.sif bcftools --version`. That output is exactly what's
needed to fix the expression.

Also run it outside the container for comparison if a system bcftools exists — a pass
outside and a fail inside pinpoints the image.

---

## 4. Launch

```bash
nextflow pull cgroza/GraffiTE      # REQUIRED: refreshes the cached asset so -r sees the new branch

nextflow run cgroza/GraffiTE -r v1.1dev-human-filter -latest \
  -profile cluster \
  --RM_dir    <RM_DIR> \
  --reference <REFERENCE.fa> \
  --human \
  --genotype false \
  --out       human_filter_test_out \
  -work-dir   <SCRATCH>/work
```

- **Use a fresh `--out`.** `publishDir` mode is `copy`, so pointing at the previous run's
  output directory leaves a stale `pangenome.trusted.human.vcf` sitting next to the new
  `pangenome.human.vcf` and makes the comparison in §5g ambiguous.
- **TSD search re-runs** and is the bulk of the runtime — `concat_repeatmask` consumes
  `tsd_report`'s output, so it cannot be skipped via `--RM_dir`. It is batched at
  `--tsd_batch_size` variants per task (default 100).
- **Faster alternative:** if the original run's `work/` directory still exists, add
  `-resume` and point `-work-dir` at it. Only `concat_repeatmask`'s script changed, so
  everything upstream — including TSD — comes from cache and the run takes minutes.
  If `-resume` re-executes TSD anyway, the cache is stale; just let the normal run
  proceed rather than fighting it.

---

## 5. Post-run checks

Work in `human_filter_test_out/3_TSD_search/`. Each check has a pass criterion — record
the actual numbers, not just pass/fail.

**a. Expected files exist, trusted files do not**

```bash
ls -la human_filter_test_out/3_TSD_search/
```
Must be present: `pangenome.vcf`, `pangenome.human.vcf`,
`pangenome.presence-absence.tsv`, `pangenome.presence-absence_human.tsv`,
`human_filter_summary.txt`, `hervk_polymorphism_summary.md`.
Must be **absent**: `pangenome.trusted.vcf`, `pangenome.presence-absence_trusted.tsv`,
`pangenome.trusted.human.vcf`.

**b. Version stamp**

```bash
grep -m1 GraffiTE_version pangenome.human.vcf
```

**c. Subfamily purity**

```bash
bcftools query -f '%INFO/matching_classes\t%INFO/repeat_ids\n' pangenome.human.vcf \
  | sort | uniq -c | sort -rn
```
Every `repeat_ids` value must start with `AluY`, `L1HS`, `SVA_D/E/F`, `HERVK-int`,
`LTR5_Hs`, `LTR5A` or `LTR5B` (possibly with an `(x)` or `(VNTR_only)` suffix). Assert it
as an allow-list — a deny-list of known-bad names misses anything unanticipated:

```bash
bcftools query -i 'n_hits==1' -f '%INFO/repeat_ids\n' pangenome.human.vcf \
  | sed 's/(.*//' \
  | grep -Evc '^(AluY|L1HS|SVA_[DEF]|HERVK-int|LTR5_Hs|LTR5A|LTR5B)'
```
Expected: `0` (any non-zero count is the number of offending records — drop the `c` from
`-Evc` to see them). Restricted to `n_hits==1` because the HERVK+SVA pair records
legitimately carry a second, non-whitelisted `SVA_A` id; those are checked in **d**.

Sanity-check that the assertion actually detects violations by running it against
`pangenome.vcf` in the same directory — that must return a large non-zero count.

**d. The HERVK+SVA carve-out fired**

```bash
bcftools view -H -i 'n_hits==2' pangenome.human.vcf | wc -l
bcftools query -i 'n_hits==2' -f '%INFO/repeat_ids\t%INFO/match_lengths\t%INFO/SVLEN\t%INFO/HERVK_CLASS\t%INFO/HERVK_PMAP\n' pangenome.human.vcf
```
Every `n_hits==2` record must be `HERVK-int` + SVA — no other multi-hit records may
appear. On the CaG cohort, 14 records are admitted by the filter and **12 survive**
`hervk_classify --strict` (the other 2 classify as `truncated_prov` at pmap ≈0.85, below
the 0.90 threshold). A different cohort will have different counts; report what you see.

**e. HERV-K annotations**

```bash
bcftools query -f '%INFO/HERVK_CLASS\t%INFO/HERVK_PMAP\n' pangenome.human.vcf \
  | grep -v '^\.' | sort | uniq -c
```
No `other`, and every `HERVK_PMAP >= 0.90` (that's what `--strict` guarantees).

**f. The polyA rule is now genuinely enforced**

The old trusted filter used `matching_classes!~"LINE"`, which does not negate on these
`Number=.` fields and so never fired. The human filter states the rule positively, so it
now actually applies:

```bash
bcftools view -H -i '(matching_classes="SINE/Alu" | matching_classes="LINE/L1" | matching_classes="Retroposon/SVA") & polyA!="TRUE"' pangenome.human.vcf | wc -l
```
Expected: `0`. `LTR/ERVK` and `Simple_repeat` records are exempt by design and may have
`polyA=FALSE` or `NA`.

**g. Before/after comparison** (only if the previous human VCF was found in §2)

```bash
OLD=<prev>/3_TSD_search/pangenome.trusted.human.vcf
NEW=human_filter_test_out/3_TSD_search/pangenome.human.vcf
bcftools query -f '%ID\n' "$OLD" | sort > /tmp/old.ids
bcftools query -f '%ID\n' "$NEW" | sort > /tmp/new.ids
comm -23 /tmp/old.ids /tmp/new.ids | wc -l    # dropped by the new filter
comm -13 /tmp/old.ids /tmp/new.ids | wc -l    # newly admitted
```
Then check *what* was dropped and gained:

```bash
comm -23 /tmp/old.ids /tmp/new.ids > /tmp/dropped.ids
bcftools query -f '%ID\t%INFO/matching_classes\t%INFO/repeat_ids\n' "$OLD" \
  | grep -Ff /tmp/dropped.ids | cut -f2,3 | sed 's/(x)//' | sort | uniq -c | sort -rn
```
The dropped set must be entirely old subfamilies; the gained set must be entirely
`n_hits==2` HERVK+SVA records. **If anything else appears in either set, that is a
finding — report it.**

Anchor numbers, if and only if the `RM_dir` is the 20-sample CaG cohort: old 5855 →
5747 after the subfamily filter → +12 surviving pairs → **5759**, with drops of 52 Alu,
32 L1, 3 SVA, 13 SVA-VNTR, 8 non-HML-2 ERVK. For any other dataset, report the observed
relationship instead of trying to match these.

**h. The filter's own report**

```bash
cat human_filter_summary.txt
```
It contains the exact bcftools expression used, the kept tally by class and subfamily, and
the dropped pME-class tally. Include the expression line and both tallies in your report.

---

## 6. Triage

| symptom | likely cause |
|---|---|
| `pangenome.trusted.vcf` still produced | Running old cached code. `nextflow pull cgroza/GraffiTE`, confirm `-r v1.1dev-human-filter -latest`. |
| 0 or absurdly few records kept | Container bcftools regex semantics — §3 should have caught this. Re-run §3 and report. |
| `no such tag defined in the VCF header` for `ULTRA_TR_span` / `polyA` | The `RM_dir` predates those INFO fields. Needs a newer RepeatMasker stage; not a filter bug. Report the `RM_dir` provenance. |
| `checkIfExists` error on `genotypes_repmasked_filtered.vcf` | `--RM_dir` points at the wrong level. It must be the `2_Repeat_Filtering` directory itself, not a numbered subdir and not the run root. |
| `n_hits==2` count is 0 | Check the cohort actually has such records before calling it a regression: `bcftools view -H -i 'n_hits==2 & matching_classes="LTR/ERVK" & matching_classes="Retroposon/SVA" & repeat_ids~"^HERVK-int"' pangenome.vcf \| wc -l` |
| TSD search takes far longer than expected | Expected — it re-runs in full. Use `-resume` against the original `work/` if it exists. |

---

## 7. Report back

Send the user:

1. Pass/fail for §3, with the container's `bcftools --version`.
2. The full `human_filter_summary.txt`.
3. Counts from §5c (purity assertion), §5d (pair count, admitted vs surviving), §5g
   (dropped / gained, with the breakdown).
4. Anything that failed, with the command and its output.

Do **not** merge this branch — it is a test branch. Do not "fix" unexpected results by
editing the filter expression in place; report them so the change can be corrected at the
source.
