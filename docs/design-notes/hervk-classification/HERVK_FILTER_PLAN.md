---
title: HERV-K classification — implementation plan
description: Design note — IMPLEMENTED — kept for the reasoning
---


!!! success "IMPLEMENTED — kept for the reasoning"
    Unlike the other design notes, this one **shipped**. It is retained because it records the reasoning behind the model. The implementation lives in `bin/hervk_classify.py`; for current behaviour see [Human mobile element insertions](../../guides/human-mei.md) and [HERV-K biology](../../background/hervk-hml2.md). Note the repository layout suggested in this plan was not followed, and the `hervk_classify_v3.py` alongside it is a superseded reference copy.

# HERV-K Insertion Polymorphism Filter — Implementation Plan

**Project**: Probabilistic post-annotation filter for HERV-K (HML-2) SV polymorphism in pangenome data
**Pipeline context**: GraffiTE v1.1dev → annotated VCF/TSV → this filter
**Activation**: Only when GraffiTE is invoked with `--human`. The filter must not run on non-human pipelines.
**Reference module**: `hervk_classify_v3.py` (in this directory) — drop-in importable, no rewriting needed.
**Regression fixture**: `HERVK_annot_test.tsv` (23 SVs) — expected outputs documented in §7.

---

## 1. Goal

Classify each TE-annotated structural variant in a pangenome VCF into one of five states based on its size and HERV-K family content:

| Label | Interpretation | Canonical \|SVLEN\| |
|---|---|---|
| `H_C` | null ↔ solo-LTR | 968 bp |
| `H_T` | null ↔ truncated proviral (single LTR + partial INT) | 1500–8000 bp |
| `H_B` | solo-LTR ↔ proviral | 8504 bp |
| `H_A` | null ↔ proviral | 9472 bp |
| `H_X` | non-transposition / SV containing HERV-K fragment | — |

Output a posterior probability per hypothesis per SV, a MAP class with its posterior, and a downstream interpretation of the per-haplotype 0/1 calls in terms of allelic states (null / solo / truncated_prov / prov). Emit these annotations into the trusted and human-filtered VCF/TSV outputs that GraffiTE produces.

A single SV record encodes a *bi-allelic* site. The 3-allele case (null/solo/proviral all segregating at one locus) appears as two overlapping SVs at the same locus — see §9.2.

## 2. Reference architecture

```
LTR5_Hs           = 968 bp     (HML-2 LTR)
HERVK-int         = 7536 bp    (HML-2 internal sequence)
solo-LTR          = 1 × LTR
proviral element  = LTR–INT–LTR  =  2 × LTR + INT  =  9472 bp
```

The five hypotheses correspond to:
- **H_C**: deletion or insertion of a single LTR (the dominant polymorphic state in humans)
- **H_T**: insertion of a 5′-truncated proviral element (one LTR + variable amount of INT, typically a 5′ deletion that still leaves the 3′ LTR intact)
- **H_B**: the LTR–INT block representing the difference between solo and proviral (LTR-LTR recombination removes 1 LTR + 1 INT = 8504 bp, converting proviral to solo)
- **H_A**: a complete proviral element (rare in humans)
- **H_X**: catch-all for SVs that contain HERV-K fragments without being insertional polymorphisms (e.g. larger SVs that incidentally span an LTR; segmental duplications carrying HERV-K material)

## 3. Mathematical model

For each SV `i`, observe `(s, λ, ν)`:

- `s` = `|SVLEN|`
- `λ` = bp matching any **LTR family** entry: `LTR5_Hs`, `LTR5A`, `LTR5B`
- `ν` = bp matching the **INT family** entry: `HERVK-int` (only)

Other LTR/ERVK family members are not pooled — they occasionally appear in RepeatMasker output but rarely correspond to true HML-2 insertion polymorphisms. SVs whose hits are exclusively non-HML-2 will fall to `H_X`, which is the correct outcome.

### 3.1 Gaussian hypotheses (H_C, H_B, H_A)

Per-hypothesis log-likelihood (up to additive constant):

```
log P(s, λ, ν | H_k) = -½ [ ((s − s*_k) / σ_{s,k})²
                          + ((λ − λ*_k) / σ_λ)²
                          + ((ν − ν*_k) / σ_ν)²
                          + ((s − (λ+ν)) / σ_t)² ]
```

Expected values per hypothesis:

| | s* | λ* | ν* |
|---|---|---|---|
| H_C | 968 | 968 | 0 |
| H_B | 8504 | 968 | 7536 |
| H_A | 9472 | 1936 | 7536 |

The fourth term `((s − (λ+ν)) / σ_t)²` enforces SV ≈ HERV-K coverage and is what penalises SVs that *contain* a fragment but aren't dominated by HERV-K sequence.

### 3.2 Truncated proviral (H_T)

A flat density on `s` over the plausibility window, with Gaussian terms on the LTR coverage (one full LTR expected) and on overall HERV-K coverage:

```
log P(s, λ, ν | H_T) = -log(T_max − T_min)               if T_min ≤ s ≤ T_max
                       − ½ ((λ − 968) / σ_λ)²
                       − ½ ((s − (λ+ν)) / σ_t)²

                     = −∞                                 otherwise
```

with `T_min = 1500, T_max = 8000`. The flat-`s` term reflects that truncation can produce any size in this range, with no preferred value. ν is left unconstrained beyond the coverage term — a truncated proviral can have anywhere from a sliver to nearly all of the INT sequence.

### 3.3 Background (H_X)

Flat over the plausible feature space:

```
log P(s, λ, ν | H_X) = -log(S · Λ · N)     S=30000, Λ=5000, N=30000
```

This sets a floor: if no structured hypothesis explains the data better than uniform, `H_X` wins.

### 3.4 Posterior

```
P(H_k | x) = π_k · P(x | H_k) / Σ_j π_j · P(x | H_j)
```

## 4. Default σ values

| Param | Value (bp) | Rationale |
|---|---|---|
| σ_{s,C} | 30 | Solo LTRs are very length-uniform |
| σ_{s,B}, σ_{s,A} | 800 | Allow ~10% INT truncation |
| σ_λ | 300 | Tight LTR annotation, but absorbs the `(x)`-merge artifact (proviral SVs often annotate as `HERVK-int(x)` with λ=0 instead of λ≈968) |
| σ_ν | 1000 | INT often partially truncated; absorbs `(x)`-merge |
| σ_t | 200 | SV must be mostly HERV-K |

## 5. Default priors

Priors are subjective and the literature on HERV-K polymorphism has known sampling biases. The recommended starting values are based on:

- **Wildschutte et al. 2016** (PNAS) catalogued ~36 polymorphic HML-2 insertions in 1000 Genomes data. Their detection method (discordant-read mapping on Illumina) is biased toward longer events, so their solo:proviral ratio overstates proviral abundance relative to what assembly-based pangenomes recover.
- **Reference-genome surveys** (e.g. Subramanian et al. 2011) consistently find solo-LTRs outnumbering full-length proviruses by ~10:1 or more in fixed insertions, reflecting the steady action of LTR–LTR recombination over evolutionary time.
- **HPRCv2 inspection by the user** shows that polymorphic full-length proviral insertions are very rare in the assembly-based pangenome, and solo-LTRs dominate the polymorphic catalog. Truncated proviral elements are more common than full-length proviral but still much less common than solo-LTRs.
- **Substantial H_X mass is needed** because pangenome SV calling produces many SVs that incidentally overlap HERV-K family sequence (large SVs spanning genomic LTRs, segmental duplications, etc.) without being HML-2 insertion polymorphisms.

Recommended defaults — applied per-SV among candidates that already have at least some HERV-K bp:

| Hypothesis | Prior |
|---|---|
| H_C (null ↔ solo) | 0.55 |
| H_T (truncated proviral) | 0.08 |
| H_B (solo ↔ proviral) | 0.05 |
| H_A (null ↔ proviral) | 0.02 |
| H_X (other) | 0.30 |

These should be exposed as configurable parameters (e.g. CLI flags or a config file), since users may want to adjust them for non-1000G/HPRC populations or for studies focused on a particular allelic state.
> I suggest to expose it as a config file to not clog the CLI `HERVK.config` or something like that 

## 6. Inputs

```
inputs/
└── <sample>.GraffiTE.annotated.tsv     # GraffiTE per-SV TSV
                                        # (only when GraffiTE was run with --human)
```

Expected GraffiTE TSV columns:
`CHROM, POS, END, ID, SVTYPE, SVLEN, n_hits, match_lengths, repeat_ids, matching_classes, fragmts, total_match_length, total_match_span, total_repeat_span` plus per-haplotype 0/1 columns named `<sample>_hap1` / `<sample>_hap2`.

A toy regression fixture is provided at `HERVK_annot_test.tsv`.
> the toy example was provided from the "trusted" tsv which only includes n_hits = 1 -- during my exploration of the Wildschutte variant in HPRC, I noted that some cases of proviral showed annotated as n_hits = 2 with HERVK + SVA --> I suspect that SVA is because SVA contains homology to LTR5_Hs. Example from pangenome.presence-absence.tsv `chr11  118740086   NA  node|15346154node|15346453_5824 NA  -8192.0 2   8190,328    HERVK-int(x),SVA_A  LTR/ERVK,Retroposon/SVA 3,1 C,C 7289,7290   8192    0.999878    None    0   0.0 0.999878    cagag` (minus the presence-absence columns). Thus it might be useful to consider these cases only (HERVK + SVA_[A-B-C-D-E-F]) as long as the total SV length is compatible with the hypothesis.

## 7. Reference module — `hervk_classify_v3.py`

The module is complete and importable. Key entry points:

| Function | Purpose |
|---|---|
| `parse_hits(match_lengths, repeat_ids)` | Returns `(λ, ν)` from comma-separated GraffiTE annotation, handling `(x)` merge suffix, using the LTR/INT family whitelists |
| `log_likelihood_gaussian(s, λ, ν, hyp, sigmas)` | H_C / H_B / H_A likelihood |
| `log_likelihood_truncated(s, λ, ν, sigmas)` | H_T likelihood |
| `classify_row(s, λ, ν, priors, sigmas)` | Posterior dict over `{C, T, B, A, X}` for one SV |
| `classify_table(df, priors, sigmas)` | Vectorised — adds `lambda_LTR, nu_INT, P_H_C, P_H_T, P_H_B, P_H_A, P_H_X, MAP_class, MAP_posterior` to the DataFrame |
| `interpret_genotype(map_class, gt_value, svtype)` | Converts a 0/1 genotype call into an allelic state string (`null`, `solo`, `truncated_prov`, `prov`, or `?`), correctly handling INS vs DEL polarity |

Module-level constants exposed for tweaking: `EXPECTED, T_MIN, T_MAX, LTR_FAMILY, INT_FAMILY, DEFAULT_SIGMAS, DEFAULT_PRIORS, ALLELE_INTERPRETATION`.

### Expected output on the regression fixture

Running `python3 hervk_classify_v3.py HERVK_annot_test.tsv` should produce:

| n | SVLEN | repeat_ids | MAP class | Posterior |
|---|---|---|---|---|
| 17 | ~968 | LTR5_Hs | H_C | 1.000 |
| 3 | 8504, 8223, 8212 | HERVK-int(x) | H_B | 1.000 |
| 1 | 5405 | HERVK-int(x) | H_T | ~0.813 |
| 1 | 3985 | HERVK-int | H_T | 1.000 |
| 1 | 318 | HERVK-int | H_X | 1.000 |

If outputs differ, something has regressed and Claude Code should investigate before proceeding.

## 8. Tasks for Claude Code

### Task A — `--human` gating
The filter must only run when `--human` is set in the GraffiTE invocation. Locate where GraffiTE branches on `--human` (likely in the main Nextflow pipeline or post-processing wrapper) and add the filter as a downstream step on that branch only. For non-human runs, the filter is skipped entirely and the trusted/human-filtered outputs do not gain HERV-K columns.

### Task B — pre-filter to single-TE-hit HERV-K SVs
Filter the GraffiTE TSV to rows where:
- `n_hits == 1` (RepeatMasker `processRepeat` already merges proviral LTR+INT into one entry tagged `(x)`, and the user wants to apply this filter on the single-TE-hit subset to keep the per-row interpretation simple)
> perhaps `n_hits <= 2` to catch the HERVK+SVA cases I was mentioning earlier
- `matching_classes` contains `LTR/ERVK`
> or `LTR/ERVK,Retroposon/SVA` if `n_hits == 2`
- After `parse_hits()`, `(λ + ν) > 0` (i.e. some bp matches the HML-2 family whitelists)

Save as `intermediate/hervk_candidates.tsv`. Report counts: total SVs in, kept, dropped (with reason).

### Task C — classify
Run `classify_table()` over the candidates with `DEFAULT_SIGMAS` and `DEFAULT_PRIORS` (or user-supplied overrides). Save `intermediate/hervk_classified.tsv`.

### Task D — emit annotations into trusted and human-filtered outputs

GraffiTE produces (at minimum) a "trusted" VCF/TSV and a "human-filtered" VCF/TSV when `--human` is set. For each:

1. **TSV output** — append columns:
   - `HERVK_class` (one of `null_solo`, `truncated_prov`, `solo_prov`, `null_prov`, `other`)
   - `HERVK_pmap` (the MAP posterior)
   - `HERVK_lambda` (LTR bp)
   - `HERVK_nu` (INT bp)
   - For SVs not in the HERV-K candidate set: leave these columns `NA`.

2. **VCF output** — add INFO fields:
   - `INFO/HERVK_CLASS=<null_solo|truncated_prov|solo_prov|null_prov|other>`
   - `INFO/HERVK_PMAP=<float>`
   - Optional FORMAT field per sample with the interpreted allele on each haplotype (e.g. `solo|null` or `prov|solo`), using `interpret_genotype()` and accounting for INS vs DEL polarity. Skip the FORMAT field if MAP posterior < 0.90 (call it `?|?` for those).

3. **VCF header** — add the corresponding `##INFO=` and `##FORMAT=` lines describing the new fields.

The trusted and human-filtered outputs should both gain these annotations identically — the filter is informational, not exclusionary. Downstream consumers can choose to filter on `HERVK_CLASS != other` and `HERVK_PMAP >= 0.9` themselves.
> I would say that the trusted and human should by default filter for `HERVK_CLASS != other` and `HERVK_PMAP >= 0.9` while leaving the other cases in the main VCF/TSV, the goal being to be conservative in the trusted/human.

### Task E — produce a summary
Save `outputs/hervk_polymorphism_summary.md` with:
- Total candidate SVs
- Per-class counts and median \|SVLEN\| per class
- Distribution of MAP posteriors (how many sites are confident vs ambiguous)
- Count of `(x)`-merged annotations and how they split across classes
- Per-haplotype allele frequency for solo-LTR, truncated_prov, and proviral states across the cohort
- List of any sites flagged as polyallelic candidates (see §9.2)
> export this summary in the same directory than the VCF and TSV

## 9. Decision points / open questions

These are flagged for the user, not for autonomous resolution:

### 9.1 H_T `T_MIN` / `T_MAX` boundaries
Defaults are 1500 / 8000 bp. The lower bound separates H_T from H_C (a 1200 bp SV is more likely a slightly-extended solo-LTR than a heavily-truncated proviral). The upper bound separates H_T from H_B (an 8200 bp SV is more likely a slightly-truncated proviral than the canonical solo-prov difference). Both are exposed as module constants and can be tuned if the empirical distribution suggests different cutoffs.

### 9.2 Polyallelic site detection
If two HERV-K-classified SVs sit within ~100 bp of each other and one has MAP class `H_C` while the other has MAP class `H_B` (or `H_T`), the locus may carry all three alleles (null/solo/proviral) at appreciable frequency. Detect these as a post-processing step on the classified table and **flag them in the summary** but do not auto-merge — the per-haplotype calls need careful interpretation that depends on which haplotype carries which SV.

### 9.3 Genotype polarity (INS vs DEL)
`interpret_genotype()` assumes the standard convention: SVTYPE=INS means alt is the longer allele; SVTYPE=DEL means alt is the shorter allele. Verify this against the GraffiTE 1.1dev VCF spec before relying on the per-haplotype allele state output. If the convention differs, flip the swap logic in `interpret_genotype()`.

### 9.4 `(x)`-merge inflation of σ_λ
RepeatMasker's `processRepeat` collapsing of adjacent LTR+INT hits into a single `HERVK-int(x)` annotation produces λ=0, ν≈8504 instead of the canonical λ≈968, ν≈7536 for proviral SVs. The default σ_λ=300 absorbs this; if a future pipeline change preserves the LTR/INT split, σ_λ can be tightened to ~100.
> the (x) logic was implemented by me in the annotation process. It may also happen (not verified) if two int subfamilies compete for the annotation, not necessarily LTR vs INT, though this is the most likely explanation as you deduced.

### 9.5 Prior tuning
The defaults in §5 are deliberately conservative for assembly-based pangenome data. Users running on heavily-biased datasets (e.g. discovery-set enriched for proviral) should override priors via the CLI/config interface.

## 10. Repository layout suggestion

```
hervk_filter/
├── README.md                          # this file
├── hervk_classify_v3.py               # core module
├── tests/
│   └── HERVK_annot_test.tsv           # regression test fixture
└── (integrated into GraffiTE main pipeline; no standalone CLI required
   since the filter is a post-annotation step gated on --human)
```

## 11. Dependencies

```
python >= 3.10
pandas >= 2.0
numpy >= 1.24
pysam                     # for VCF I/O when emitting Task D outputs
```

No R dependency for the core; an R port of the pre-H_T classifier exists at `hervk_classify.R` for ad-hoc downstream analysis but is not part of the production filter.

# 12. Populate the README 

Document this extensively in a dedicated section of the README for HERVK

---

## Quick-start for Claude Code

```bash
# 1. Confirm the module passes regression on the toy fixture
cd hervk_filter/
python3 hervk_classify_v3.py tests/HERVK_annot_test.tsv
# -> expected output documented in §7

# 2. Locate the --human branch in GraffiTE 1.1dev and wire in:
#    - Task B (pre-filter)
#    - Task C (classify)
#    - Task D (emit into trusted + human-filtered VCF/TSV)
#    - Task E (summary)
```

The first concrete deliverable is the integration into GraffiTE's `--human` branch with the four task steps wired in. Stop after Task E and surface the summary for review before considering further refinements (prior tuning, polyallelic merging).
