---
title: HERV-K classifier v2 — evidence-first allele states
status: implemented (stage 3); stage 4 consolidation deferred to phase 2
supersedes: HERVK_FILTER_PLAN.md (kept for its reasoning; do not implement from it)
---

# HERV-K (HML-2) classifier v2

!!! note "Why this exists"
    v1 decided REF and ALT allele states from expected size arithmetic alone.
    That cannot separate an 8.5 kb LTR+INT block entering a solo LTR from a
    complete provirus entering an empty site, and it could not reach
    `null_prov` at all. This note records the evidence that forced the
    redesign, so the reasoning survives the code.

## 1. What was wrong

`bin/hervk_classify.py` v1 ran a Gaussian MAP over `(|SVLEN|, λ, ν)` against
canonical sizes 968 / 8504 / 9472, where λ is LTR5 bp and ν is `HERVK-int` bp.
It never looked at the reference genome and never saw the internal structure of
the variant allele. Two independent lines of evidence say it was systematically
wrong in one place:

1. **Reference inspection.** At `solo_prov` loci that were checked by hand,
   CHM13v2 carried no solo LTR — the reference looked *null*. v1 had asserted
   `REF = solo` purely because `|SVLEN| ≈ 8504`.
2. **Wildschutte 2016.** The two non-MEI loci that matched in the liftOver
   comparison, `chr19-21797327-INS-9478` and `chr19-22370220-INS-8229`, are
   labelled `pro_pre` — provirus versus pre-integration site — while GraffiTE
   called both `solo_prov`.

Genome-wide `null_prov` was 0, which is not plausible. Reaching it needs
λ ≈ 2 × 968 = 1936 bp, and λ is structurally never that large.

## 2. Root cause: the architecture was there all along

`bin/annotate_vcf.R:73-88` groups RepeatMasker fragments by link ID and
collapses each group to the top-SW-score family name plus `(x)`, with the span
set to the full extent. That is the right summary for most TE classes. For
HML-2 it erases exactly what matters: a full provirus arrives as
`HERVK-int(x)` covering the whole SV, with **zero** LTR bp.

The raw tables are intact. Reading
`2_Repeat_Filtering/*/repeatmasker_dir/indels.fa.out` directly recovers the
structure (query coordinates → consensus coordinates, element orientation):

```
chr19-21797327-INS-9478   LTR5_Hs[1-968] · HERVK-int[1-7536] · LTR5_Hs[1-968]
chr8-7226885-INS-9468     LTR5_Hs[1-968] · HERVK-int[1-7536] · LTR5_Hs[1-968]
    two full terminal LTRs -> a complete provirus -> REF must be null

chr6-78894876-INS-8465    LTR5_Hs[411-968] · INT · LTR5_Hs[1-410]
chr11-101705464-INS-8498  LTR5_Hs[575-968] · INT · LTR5_Hs[1-574]
chr7-4700334-INS-8504     SVA_A[951-1113]  · INT · LTR5_Hs[175-968]
chr12-58305931-DEL-8489   SVA_A[940-1242]  · INT · LTR5_Hs[304-968]
chr12-133148145-INS-4933  LTR5A[736-1033]  · INT · LTR5A[1-735]
    one LTR, split complementarily across the termini -> REF = solo
```

So `bin/hervk_arch.py` reads the raw tables. `annotate_vcf.R` is left alone —
it feeds every other class and its summary is correct for them.

## 3. The two signatures

**`ARCH_2LTR`** — both terminal LTRs cover the full consensus (~1..L). Nothing
was consumed by the alignment, so the SV carries the whole element and the
reference was empty. `REF = null`, `ALT = provirus`.

**`ARCH_PERM`** — one LTR is split across the two termini, the consensus
intervals complementary: the 5′ fragment covers `[k+1..L]` and the 3′ fragment
covers `[1..k]`, summing to one consensus length.

The mechanism: when the reference holds a solo LTR, that LTR matches *both*
LTRs of the provirus. The aligner must break it somewhere; reference positions
`1..k` anchor to the alt's 5′ LTR and `k+1..L` to the alt's 3′ LTR, so what
falls between the anchors is `LTR[k+1..L] · INT · LTR[1..k]`. The inserted
length is `(L−k) + INT + k` — constant for every `k`.

### `k` is an alignment property, not a biological one

The recombinant solo LTR's own crossover point is unrecoverable (the parental
LTRs are near-identical) and irrelevant. Because the inserted length is the
same for every `k`, the gap penalty cancels, and the placement is decided by
the handful of substitutions separating the two LTR copies from the reference
LTR. Measured in CaG:

- **`k` varies between haplotypes at one locus.** chr11 is called at `k = 0`
  from sample 11118731 and `k = 574` from sample 11107362 — both PAV, both
  `CALL_SOURCE=CIGAR`.
- **Two independent measurements agree.** Record B sits 574 bp from record A,
  and its RepeatMasker permutation point is consensus 575, i.e. `k = 574`. The
  chr1 pair repeats this: 871 bp apart, `k = 868`.
- **Left-alignment is not the driver.** `bcftools norm -f` in `truvari_merge`
  moved **1 of 52** HERV-K records; the other 51 shifted only by the 1 bp
  anchor-base convention. minimap2's DP has usually already parked the
  breakpoint against a mismatch.
- **Observed `k`:** 0 (×7), then 102, 174, 303, 410, 574, 735.

Nothing downstream may key on `k`, and it must not be read as a measurement of
the recombination event. `k` is also expected to be caller-dependent — this is
**untested**, since every CaG record is PAV.

## 4. The `SVA_A` hit is not a mimic

v1 treated an SVA hit of 250–400 bp beside HML-2 as "mimicking" LTR5_Hs. The
real relationship is homology: SVA's SINE-R domain is HERV-K LTR derived. Every
SVA hit seen beside an HML-2 element in CaG sits at SVA consensus ≥ ~900, and
where it is the *sole* annotation of a terminus its length completes the LTR
exactly — chr12-58305931 is 303 bp of SVA plus `LTR5_Hs[304-968]`, summing to
968. v1's window got the right answer for the wrong reason and missed
chr7-4700334 (172 bp) and chr11-101705464 (61 bp).

v2 reassigns SVA hits at consensus ≥ 900 abutting HML-2 to the LTR class, and
tiles hits winner-take-all on the query axis so an SVA hit lying *on top of* an
LTR5_Hs hit is counted once rather than twice.

## 5. The reference check is the workhorse, not a fallback

Sequence alone cannot resolve the degenerate case: one full LTR at one terminus
and nothing at the other. That is either a permutation at `k = 0` (insertion
right at the LTR boundary, `REF = solo`) or a one-LTR-truncated provirus into
an empty site (`REF = null`). **Roughly half the CaG candidates land in that
form**, so `bin/hervk_ref_state.py` — which masks a window of the reference at
each candidate and reads off what is there — carries about half the loci.
`ARCH_PERM` is a corroborating shortcut, not the primary mechanism.

`chr19-22370220` is the case that proves it: LTR + INT with a 291 bp internal
deletion, 8229 bp total. Size says LTR+INT-into-solo; the reference says null,
and Wildschutte says `pro_pre`.

## 6. Resolution order

| Rule | Condition | REF | ALT |
|---|---|---|---|
| `ARCH_2LTR` | two full terminal LTRs | `null` | `provirus` |
| `ARCH_PERM` | one LTR, complementary split | `solo` | `provirus` |
| `ARCH_SOLO` | lone LTR, no internal region | INS `null` / DEL `solo` | INS `solo` / DEL `null` |
| `REF_ANNOT` | degenerate — the masked reference decides | from the reference | REF ± SV content |
| `DENOVO_LTR` | *(deferred, see §9)* | | |
| `UNRESOLVED` | nothing resolved it — kept and flagged | `.` | `.` |

`HERVK_PMAP` is now a properly normalised confidence over the resolved class.
It is reported and used by `--strict`; it never decides anything. (v1 mixed
unnormalised Gaussians of different σ with normalised flat densities, which
silently favoured the narrow solo-LTR hypothesis by ~26×.)

## 7. Locus grouping

Records describing the same insertion can be anchored anywhere inside the
shared reference LTR, so they land up to **one LTR apart**. That is the
derivable clustering radius — `hervk_locus_window = 1200` — not an arbitrary
window. Observed offsets in CaG: 107 bp (chr12), 559 bp (chr6), 574 bp
(chr11), 871 bp (chr1). Three of those exceed truvari's default 500 bp
`refdist`, which is why truvari did not collapse them.

Reordering `bcftools norm` before `truvari collapse` would **not** fix chr11:
collapse currently sees 823 bp, norm would bring it to 574 bp, still above 500.

Grouping runs over **every** candidate in `pangenome.vcf`, not the human
subset. `--human` requires `FILTER="PASS"` and PAV emits `TRIM`/`COMPOUND` on
real HML-2 records — `chr15-2092086-DEL-8221` is one — so grouping from the
human VCF alone could split a locus whose partner was removed for an unrelated
reason. Only the human VCF is *written*; `pangenome.vcf` induces the graph and
stays byte-identical.

## 8. Result on the CaG candidates

- `null_prov` becomes non-zero for the first time: `chr19-21797327` and
  `chr8-7226885`, both previously `solo_prov`. The first matches Wildschutte's
  `pro_pre`.
- All six predicted permutation points reproduce exactly.
- Four multi-record loci flagged: chr12:55,299,985 (`MULTIALLELIC` — the
  three-allele locus), chr6:78,894,316 (`POLARITY_CONFLICT`),
  **chr11:101,704,640** (missed by snarl-based grouping, which keys on
  `INFO/RC`/`RD` and sees two different snarls), and chr1:75,219,429.
- 80 non-HML-2 `LTR/ERVK` records (`HERVK9-int`, `MER11A`, `LTR13` …) are now
  labelled `other`. v1 skipped them silently, which also let them pass
  `--strict` unfiltered.

**chr6:78,894,316 is unresolved as a biological question.** The DEL implies
`REF = provirus`; the INS implies `REF = solo`. Both cannot hold. The DEL's REF
field is literal reference sequence containing a full LTR + full INT, so the
masked reference should say `provirus` — which would invert the solo/provirus
frequencies that `paper/LATEST_GT/build_hervk_vcf.py` currently reports for
that locus. Confirm against the reference before trusting either.

## 9. Deliberately deferred

- **`DENOVO_LTR`** (rule 5) is not implemented. It was specified as a fallback
  for when the reference is unavailable or ambiguous. Every degenerate CaG
  record is resolved by `REF_ANNOT`, so it buys nothing yet — but it is a real
  gap for cohorts without a good reference annotation, and the evidence code is
  reserved.
- **Stage 4 consolidation** (`hervk_reconcile.py consolidate`) — phase 2. Graph
  genotyping cannot be re-run for this cohort, so it is built against a
  supplied giraffe VCF.
- **graphaligner and pangenie** reconciliation — phase 3.
- **Sex-aware ploidy and PAR.** GraffiTE has no native handling; the hardcoded
  `-R chrX:1,chrY:1` in `vg call` calls chrX/chrY haploid for every sample
  regardless of sex. The locus layer flags non-autosomes
  (`PLOIDY_UNVERIFIED`) and otherwise stays out of the way.
- **Caller dependence of `k`** — untested; needs the route coverage test set.
