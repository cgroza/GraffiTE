---
title: Changelog
description: Release history of GraffiTE, from the first beta to the current v1.1 development branch.
---

# Changelog

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](getting-started/v1.0-vs-v1.1.md).

Entries before v1.1 are the README's, kept with their original dates. Since beta 0.2.5 releases
are identified by commit; give the commit you run when asking for help.

---

## v1.1 (development)

The `v1.1dev` branch. A code update and an image update are both needed to see every change.

**Annotation**

- Tandem repeats annotated with [ULTRA](https://github.com/TravisWheelerLab/ULTRA) after
  RepeatMasker: `ULTRA_TR` (non-redundant bp of tandem repeat) and `ULTRA_TR_span` (fraction of
  the variant, capped at 1).
- New filter metric `total_repeat_span`, the non-redundant union of RepeatMasker TE hits and
  ULTRA tandem repeats over the variant length, in place of `total_match_span` alone. Long
  insertions where a TE sits next to a polyA or tandem tail are kept. `--repeat_span_cutoff`
  (default `0.80`) sets the threshold.
- polyA tail detection: `polyA=TRUE/FALSE/NA` on single-hit records, after trimming an exact
  TSD copy from the variant end.
- `--mammal` discontinued. L1 5′ inversion detection and SVA VNTR-only handling run on every
  dataset.
- L1 5′ inversions reported as `INFO/L1_5PINV` (the hit ID of the inverted fragment, or `None`)
  in place of `mam_filter_1=5P_INV`.
- SVA hits that fall inside the VNTR region are renamed with a `(VNTR_only)` suffix and
  reclassified from `Retroposon/SVA` to `Simple_repeat`; `mam_filter_2` removed.
- OneCodeToFindThemAll removed from the annotation step; the RepeatMasker `.out` is read
  directly.
- `concat_repeatmask` accepts a gzip-compressed `--reference` and re-compresses it to BGZF.
- Satellite and RNA hits are filtered out of `pangenome.vcf`; the repeat-class filter in the
  annotation step was removed.

**Subsets**

- `pangenome.trusted.vcf` and matching `*.presence-absence*.tsv` tables next to
  `pangenome.vcf`; trusted records carry `FILTER=PASS`. Parameters `--trusted_min_svlen`
  (250), `--trusted_max_ultra_span` (0.6), `--trusted_ignore_filter`.
- `--human` rebuilt from `pangenome.vcf` by its own filter, output `pangenome.human.vcf` (was
  `pangenome.trusted.human.vcf`); `pangenome.trusted.vcf` is not written on `--human` runs. The
  filter keeps `AluY*`, `L1HS`, `SVA_D/E/F` and the HML-2 lineage, and admits the
  `HERVK-int`+`SVA` proviral pattern up to `--hervk_pair_max_hits` hits (PR #98). Whitelists
  and thresholds are parameters; `human_filter_summary.txt` reports what was kept and dropped.

**HERV-K (HML-2) classifier, `--human` only** (PRs #94 to #97)

- Evidence-first allele-state classifier from the element's architecture and from masking the
  reference at the locus (`hervk_classify.py`, `hervk_arch.py`, `hervk_ref_state.py`), with
  the `ARCH_INT_PERM` internal-region permutation signature and copy-number allele states.
- Loci clustered by overlap and flagged (`HERVK_LOCUS`, `HERVK_MEI`, `HERVK_CNV`, ...); graph
  genotypes withheld per allele at copy-number and tandem-duplication loci
  (`--hervk_mask_graph_gt_at_cnv`).
- Consolidated VCFs, one record per locus: `pangenome.human.consolidated.vcf` in discovery and
  `GraffiTE.merged.genotypes.human.vcf.gz` after genotyping, with reports; `--hervk_reconcile`,
  `--hervk_reconcile_vcf` to consolidate an existing genotyped VCF with `--genotype false`.
- `hervk_candidates.vcf` publishes the full candidate set; `--hervk_max_svlen`,
  `--hervk_ref_flank`, `--hervk_locus_window`, `--hervk_ref_annotation`, `--hervk_strict`,
  `--hervk_config`.
- `human_hervk_ids` accepts both `HERVK-int` and Dfam's bare `HERVK` for the internal region.

**Inputs and genotyping**

- `--svs`: per-sample VCFs fed straight into the merge.
- `--pav`: the PAV caller as a discovery entry point.
- `vg call` receives `-a -A` on every graph method, and ploidy 1 on `chrX` and `chrY`.
- Precomputed graph alignments: `--graph_alignments` (CSV `sample,gaf,pack`), published to
  `GraffiTE_alignments/`.

**TSD chain**

- `prepTSD.sh` reads the flanks through htslib, re-compresses a gzip reference to BGZF, and
  stops the run instead of writing an empty flank file; a missing `exact_match.py` is an error
  rather than no hits; `add_polyA.py` trims the two-copy `TSD` value correctly. Test:
  `test/tsd/test_tsd_chain.sh` (commit `884afa8`).

**Fixes on the documentation branch**

- `--tsd_win` reaches the matcher; scores were computed against a fixed 30 bp
  (`4c8e385`).
- PanGenie genotypes are split to one ALT per record before the merge; the normalized file
  was written but never emitted (`46c0726`).
- `make_graph`'s graphaligner branch interpolated `$PWD` as a Groovy variable and could not
  run (`9b77788`).
- `--vcf` beside a discovery flag, `--graph_method precomputed` without its inputs, and
  `--graffite_vcf --human` are refused with a message instead of failing on an undefined
  channel (`159de3f`, `b6a9487`, `9bc505b`); an empty `panmethyl/` submodule is reported
  (`959044b`).
- `svs`, `graph`, `graph_alignments`, `vcfs`, `lifted` and `bed` are declared in
  `nextflow.config` (`4187750`); `hervk_annotate` and `hervk_reconcile` have resource
  parameters (`f0698ec`); `manifest.version` matches `version.txt` (`e9f4f74`);
  `merge_VCFs` no longer passes `publishDir` an unknown option (`4b917b4`).
- `GraffiTE.def` installs truvari, ULTRA, pyfaidx and pypy3, which the pipeline calls and the
  published image already had (`9cfaf15`, recipe not rebuilt).
- Six unused scripts removed from `bin/` (`ebaf837`).

---

## 01/01/25

Commit [1cbebbf](https://github.com/cgroza/GraffiTE/commit/1cbebbfc0f4ccc5436670d9aa8d2023a90f2eeef).

- Removed `--nolow` from the RepeatMasker call. It could produce spurious hits on
  low-complexity regions of some TE consensus sequences, mistaking tandem repeats for TEs. Do
  not use `--nolow` with RepeatMasker outside debugging.

## 11/07/24

Commit [76537f9](https://github.com/cgroza/GraffiTE/commit/76537f9b5da4024ba03f760f58b024a5f485bf7a).

- New `--tsd_time` option for the time request of the TSD processes under the `cluster`
  profile. Default `1h`. Code update only, no new image.

## 10/22/24

Commit [47ad044](https://github.com/cgroza/GraffiTE/commit/47ad04469e475e9dcbfd4ffc17faa4ba42c5d94d).
Pull request by [Han-Cao](https://github.com/Han-Cao).

- Faster annotation of large VCFs.
- Coordinates in the SVA-VNTR module changed from 1-based to 0-based. Code update only.

## 10/21/24

- RepeatMasker coordinates transformed from 1-based to 0-based to match the BED standard and
  measure hit lengths exactly. Fixes [issue #43](https://github.com/cgroza/GraffiTE/issues/43).

## 06/24/24

- New option `--break_scaffolds`, which splits contigs at runs of more than four N. With some
  scaffolded genomes minimap2 fails with `[E::parse_cigar] CIGAR length too long at position
  ...`, a limit of htslib and the SAM specification; breaking scaffolds at N stretches avoids
  it.

## 06/17/24

- New compatible class names in the TE library: `MITE`, `TIR` and `IS`, as in `>TEnameX#MITE`,
  `>TEnameY#TIR/Mariner` or `>TEnameX#IS`. Earlier versions discarded them.
- Accepted `Class` names in `>TEname#Class/Superfamily`: `LINE`, `LTR`, `SINE`, `RC/Helitron`
  (treated as `DNA/RC`), `DNA`, `TIR`, `MITE`, `Retroposon`, `IS`, `Unknown`, `Unspecified`.
  An entry without a class is treated as `Unknown`. Any name and superfamily is accepted when
  the class is one of these.

## 02/13/24

- Versions are identified by commit ID from here on.
- The L1 inversion flag (`--mammal`) had stopped working; fixed.
- Winnowmap available as an alternative to minimap2 with `--aligner winnowmap`.

## beta 0.2.5 (09-11-23)

- Fixed a VCF annotation error when two distinct variants shared the same `POS`; annotations
  are now distinct by variant sequence.
- Cleaner GraphAligner VCF outputs.

## beta 0.2.4 (06-27-23)

- Refactored to Nextflow DSL2.

## beta 0.2.3 (02-21-22)

- SV discovery from assemblies and long reads together; the calls from both are merged before
  filtering and genotyping.
- Parameters with defaults to control time, CPU and memory for each process, for
  `-profile cluster`.
- Variants are merged only within the same `SVTYPE`.

## beta 0.2.2 (02-01-22)

- `sniffles2` as an alternative to `svim-asm`, from long reads (`--longreads`) instead of an
  assembly. At this release the two callers could not be combined.
- `--asm_divergence <asm5/asm10/asm20>` sets the minimap2 preset ahead of `svim-asm`; default
  `asm5`.
- `time`, `cpu` and `memory` directives per process.

## beta 0.2.1 (11-30-22)

- `--RM_vcf` and `--RM_dir` to start at the TSD search from the output of a RepeatMasker run
  (`2_Repeat_Filtering`), for runs that crashed in the TSD search and cannot be resumed.
- TSD search in batches of 100 variants, which divides the number of temporary directories by
  100 and spares the inode quota; batches run in parallel.

## beta 0.2 (11-11-22)

- Two new read aligners for genotyping, [giraffe](https://github.com/vgteam/vg#mapping) and
  [GraphAligner](https://github.com/maickrau/GraphAligner): `--graph_method
  [pangenie/giraffe/graphaligner]`, default `pangenie`.
- `--vcf`: a sequence-resolved VCF as input, skipping genome alignment.
- `--graffite_vcf`: a `3_TSD_search/pangenome.vcf` from an earlier run, skipping everything
  but read mapping.
- Dropped the `biomartr` dependency.

## beta 0.1 (11-02-22)

- First release.
