---
title: CLI ergonomics and benchmark proposal
description: Design note — PROPOSED — not implemented
---

!!! danger "PROPOSED — not implemented"
    The proposals in sections 2 and 3 are **not implemented**: there is no `--input` samplesheet, no `--realign`, no `pangenome.discovery_genotypes.vcf` and no `test/benchmark/` harness. Section 1, which maps the pipeline into three stages, *is* an accurate description of the code and is the basis for [Choosing your inputs](../getting-started/choosing-your-inputs.md).

# GraffiTE — CLI ergonomics & benchmarking proposal

> **Status:** Draft for colleague feedback — *not yet implemented*.
> **Authors:** Clément Goubert (requests) + Claude (analysis/drafting).
> **Date:** 2026-06-26 · **Branch:** `v1.1dev`
> **Scope note:** Excludes the `panmethyl` / epigenome (`--epigenomes`) capabilities by request.

---

## 0. What Clément asked for (summary of requests)

GraffiTE has grown into many combinations of tools depending on the input data, which is
getting confusing for casual users. Goal: **improve CLI ergonomy, command names, and sensible
defaults to reduce confusion, while keeping full customization available.** Concretely:

1. **Map all the pipelines that make sense** (the valid combinations of discovery × annotation × genotyping).
2. **Suggest CLI / quality-of-life improvements** that make the CLI simpler and more intuitive. Specifically called out:
   - Allow **extracting reads from a BAM and re-aligning** before SV search.
   - **Report the genotypes produced by sniffles2 and PAV**, collated with the TE annotation
     (these callers emit genotypes, but GraffiTE currently does *not* surface them next to the
     annotation because they are not graph-genotyped). Users should be able to compare them with
     graph genotypes, use them directly if they want to **skip graph alignment**, or **recover SVs
     that can't / fail to be graph-genotyped**.
   - Balance new features against keeping the CLI simple and intuitive.
3. **Propose a focused benchmark** (not the full combinatorial explosion) selecting a few
   detection + genotyping modes to compare against competing software. Available benchmark data:
   **HG002 / hg38**, with **HiFi reads, an HG002 assembly, and a truth VCF**. The benchmark must
   evaluate (1) **TE detection** (presence/absence) and (2) **TE genotypes**.

---

## 1. Map of all sensible pipelines

GraffiTE is **three stages in series**. The confusion comes from input flags encoding
*{data type × tool}* instead of *{what the user has}*.

### Stage A — SV discovery → `1_SV_search/SVs.vcf`
All backends produce per-caller VCFs that are mixed and merged by `truvari_merge`.
**A1–A5 can be freely combined** (e.g. assemblies + long reads together).

| Flag | Input you have | Tool chain | minimap2? | Code |
|------|---------------|------------|-----------|------|
| `--longreads` | long-read **FASTQ** (+ `type`) | minimap2/winnowmap → sniffles2 | ✅ | `main.nf:56-60`, `module/main.nf:40` |
| `--bams` | long-read BAM **already aligned to `--reference`** | sniffles2 only | ❌ | `main.nf:62-65` |
| `--assemblies` | contigs / haplotype FASTA | minimap2 → svim-asm | ✅ | `main.nf:74-82`, `module/main.nf:140` |
| `--pav` | assembly haplotypes | PAV | ✅ (internal) | `main.nf:84-88`, `module/main.nf:106` |
| `--svs` | per-sample SV VCFs (CSV) | passthrough → truvari | — | `main.nf:90-93` |
| `--vcf` | one sequence-resolved merged VCF | bypasses discovery | — | `main.nf:111` |

### Stage B — Annotation → `pangenome.vcf` (+ trusted / human / presence-absence TSVs)
RepeatMasker + ULTRA + TSD search. **Always runs**, except:
- `--graffite_vcf` → skip discovery **and** annotation (jump straight to genotyping). `main.nf:130`
- `--RM_dir` → resume from a previous RepeatMasker run. `main.nf:102`

### Stage C — Genotyping → `4_Genotyping/GraffiTE.merged.genotypes.vcf`
Controlled by `--genotype` (default `true`) + `--genotype_with` CSV + `--graph_method`:

| `--graph_method` | Read type | Engine |
|------------------|-----------|--------|
| `pangenie` (default) | short | PanGenie (k-mer) |
| `giraffe` | short **or** long | vg giraffe → vg call |
| `graphaligner` | long | GraphAligner → vg call |
| `precomputed` + `--vcfs` | — | reuse existing vg-call VCFs |
| `--genotype false` | — | stop after annotation |

### Valid space
`(one or more of A1–A5) → B → (0 or 1 of C, matched to read type)`.
`--vcf` and `--graffite_vcf` are exclusive shortcuts that bypass earlier stages.

```
            DISCOVERY (mix freely)          ANNOTATE         GENOTYPE (pick ≤1)
 fastq  ─ minimap2 ─ sniffles2 ┐
 bam    ──────────── sniffles2 ┤
 asm    ─ minimap2 ─ svim-asm  ┼─ truvari ─ RepeatMasker ─┬─ pangenie    (short)
 asm    ──────────── PAV       ┤            +ULTRA+TSD    ├─ giraffe      (short/long)
 svs    ────────────────────── ┘            = pangenome  ├─ graphaligner (long)
 vcf    ─── bypass ────────────────────────►            └─ none
 graffite_vcf ─── bypass discovery+annotate ──────────────►
```

---

## 2. CLI / quality-of-life improvements

Ordered by impact-to-effort. Intended grouping: a small "ergonomics" PR for (a,b),
then feature work for (e), then the unified samplesheet (c,d), then polish (f,g).

### High impact, low effort

**a. Kill silent-flag failures.** Nextflow does *not* normalize dashes to underscores. The README
shows `--genotype-with reads.csv`, but the param is `genotype_with` (`nextflow.config:36`), so that
command silently falls back to the default `reads.csv` and only "works" if the file is named that.
Same class of bug as `--bam` vs `--bams`. Fix:
- Accept singular/plural & dash/underscore aliases (`--bam`/`--bams`, `--longread(s)`,
  `--assembly`/`--assemblies`, `--genotype-with`/`--genotype_with`).
- Add a **strict unknown-param check** at workflow start (error + "did you mean…?"). This alone
  would have caught `--bam` immediately instead of the cryptic downstream error.

**b. Replace the generic discovery error** (`main.nf:114`, `"No --longreads, --assemblies…"`) with a
guided message: short cheatsheet of input modes + closest match to what was passed.

### High impact, medium effort

**c. One discovery input, autodetected — `--input samplesheet.csv`.** Collapse
`--longreads/--bams/--assemblies/--pav/--svs` into a single sheet with a `kind` column
(`reads|bam|asm|sv_vcf`); route each row to the right backend. Keep old flags as thin aliases for
back-compat. Biggest clarity win — users describe *what they have*, not which tool to run.

```csv
sample,kind,path,type
HG002,reads,HG002_HiFi.fq.gz,hifi
HG002,asm,HG002.hap1.fa,
```

**d. BAM → realign (requested).** Add `--realign` (or `kind: ubam`) so a BAM (unaligned HiFi uBAM, or
aligned to a *different* reference) is `samtools fastq`-extracted and re-mapped with minimap2 before
sniffles2. `bam_to_fastq` already exists in the genotyping path (`module/main.nf:485`) — promote it
into discovery. Decision rule: `--bams` = trust existing alignment; `--realign` = strip & remap.

**e. Surface discovery genotypes collated with TE annotation (requested).** Today `truvari_merge`
does `bcftools annotate -x INFO` + `+setGT -t . -n 0` (`module/main.nf:195,202`), flattening caller
GTs to 0/0 so they are never reported next to the annotation. Proposal:
- Carry the **original per-sample GTs** from sniffles2 / svim-asm / PAV through annotation into a new
  output `pangenome.discovery_genotypes.vcf` (+ presence-absence TSV), distinct from the graph-based
  `GraffiTE.merged.genotypes.vcf`.
- Enables: (i) skip graph genotyping entirely, (ii) compare caller-GT vs graph-GT, (iii) **recover
  variants that fail graph genotyping** (e.g. dropped by `vg call`).
- Keep a `GT_origin`/`SOURCE` tag for provenance (sniffles2 GTs carry quality; PAV gives diploid
  phased GTs). Gate with `--report_discovery_genotypes` (or auto-emit when GT-producing callers run).

### Medium impact

**f. `--mode` presets** wrapping the README's GT-sn-GA-style matrix, e.g. `--mode longread_hifi`
(= `graph_method=graphaligner` + sniffles defaults), `--mode shortread` (= pangenie). Presets set
defaults; every underlying flag stays overridable — one knob for casual users, full control retained.

**g. `--check` dry-validate path** that parses sheets, verifies files exist, confirms read `type` ↔
`graph_method` compatibility (e.g. error early on `pangenie` + long reads), and prints the resolved
plan — before submitting any Slurm job.

### Suggested sequencing
(a,b) guardrails → (e) discovery genotypes → (c,d) unified sheet + realign → (f,g) polish.

---

## 3. Benchmark plan (HG002 / hg38; HiFi + assembly + truth VCF)

Pick a few representative modes per question, not the full cross-product. Two axes:
**detection** and **genotyping**.

### Prep the truth set
Truth VCF is not TE-specific → restrict it: keep INS/DEL ≥ ~100 bp, intersect with a TE annotation
(RepeatMasker on truth ALT alleles, or lift the DFAM library) → `truth.TE.vcf`, stratified by family
(Alu / L1 / SVA). All comparisons use `truvari bench` (sequence-aware).

### Axis 1 — Detection (precision / recall / F1, per family)
| Mode | GraffiTE call | Tests |
|------|--------------|-------|
| `asm-svim` | `--assemblies HG002.fa` | assembly-based detection |
| `asm-pav` | `--pav` | PAV vs svim-asm |
| `reads-sniffles` | `--longreads HiFi` | read-based detection |
| `asm+reads` | both | does merging improve recall? |

**Competitors (HiFi/assembly-capable):** TLDR, PALMER (long-read MEI callers); plus *raw* sniffles2 /
svim-asm with **no TE filtering** (raw-vs-GraffiTE delta = TE-specificity gain). MELT/xTea only if a
short-read set is added later.

### Axis 2 — Genotyping (GT concordance vs truth)
`truvari bench` then genotype concordance on the matched set:
| Mode | GraffiTE call | Tests |
|------|--------------|-------|
| `GT-graphaligner` | discovery → `--graph_method graphaligner` + HiFi | graph GT, long reads |
| `GT-giraffe` | `--graph_method giraffe` + HiFi | giraffe vs GA on long reads |
| `GT-discovery` | **feature (e)**, sniffles2 direct GT, no graph | direct vs graph GT |
| `GT-pangenie` | only if a short-read set is added | short-read baseline |

Report: GT concordance, non-ref concordance, and **fraction of detected variants left ungenotyped**
by each engine (motivates feature **e**).

### Harness
Small reproducible sub-workflow in `test/benchmark/`: modes CSV → run GraffiTE → `truvari bench` →
collate `benchmark_summary.tsv` (mode, family, P/R/F1, GT-concordance). Doubles as chr-subset CI.

---

## 4. Pick-up notes (for resuming implementation)

**Decision still open:** which piece to build first. Candidates, in suggested order:
1. **(a,b) CLI guardrails** — param aliases + strict unknown-flag check + guided error. Smallest, prevents the `--bam` class of bug. Touch: top of `main.nf`.
2. **(e) Discovery genotypes** — new `pangenome.discovery_genotypes.vcf`/TSV. New capability the benchmark depends on. Touch: `truvari_merge` GT handling (`module/main.nf:158-208`), annotation passthrough, new output in `concat_repeatmask`.
3. **(c,d) Unified `--input` samplesheet + `--realign`** — biggest clarity change. Touch: input parsing block `main.nf:52-93`; reuse `bam_to_fastq` (`module/main.nf:485`).
4. **(3) Benchmark harness** — `test/benchmark/` scaffold.

**Key code anchors:**
- Input routing: `main.nf:46-95`
- Discovery-GT loss: `module/main.nf:195` (`annotate -x INFO`) + `:202` (`+setGT -t . -n 0`)
- Genotyping branch: `main.nf:135-221`
- Existing BAM→fastq: `module/main.nf:485-495`
- Final genotype merge: `merge_VCFs`, `module/main.nf:551`

**Constraints to respect:** keep old flags working as aliases (back-compat); minimap2 cannot read
BAM (drives feature d); `--graph_method pangenie` is short-read only.
