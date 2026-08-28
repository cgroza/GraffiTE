## HERV-K (HML-2) classifier v2 — evidence-first allele states

Replaces the size-arithmetic HERV-K classifier. The old model decided REF and
ALT allele states from `|SVLEN|` alone, which cannot separate an 8.5 kb
`LTR+INT` block entering a solo LTR from a complete provirus entering an empty
site — and made `null_prov` unreachable, because that class needs ~2 LTRs of
LTR5 bp and `annotate_vcf.R` collapses RepeatMasker fragments to one name plus
`(x)`, leaving a provirus with zero LTR bp.

The fix reads the architecture out of the **raw** RepeatMasker tables and, where
that is degenerate, **masks the reference** at the locus.

### What lands

| | |
|---|---|
| `bin/hervk_arch.py` | query-axis tiling, SVA SINE-R reassignment, the `ARCH_2LTR` / `ARCH_PERM` signatures |
| `bin/hervk_ref_state.py` | masks a reference window per candidate; adaptive re-cut when an element hits the window edge |
| `bin/hervk_classify.py` | rewritten around an evidence ladder; `HERVK_PMAP` is now confidence only and never decides a class |
| `bin/hervk_reconcile.py` | `flag` groups records into loci; `consolidate` resolves genotypes by allele dosage across a locus's members |
| `module/main.nf`, `main.nf`, `nextflow.config` | two new `--human`-only processes; **`pangenome.vcf` is never written to** and `vcf_ch` is not rebound, so the graph is unchanged |

### Results on CaG

- **`null_prov` is non-zero for the first time** — `chr19-21797327`, `chr8-7226885`, `chr19-22370220`.
- **12/12 concordant with Wildschutte 2016**, reached independently of it. The two `pro_pre` loci now line up under the label's plain reading; previously the comparison only agreed under an element-form reading it flagged as unverified.
- **chr11:101,704,640 is one locus, not two** — 574 bp apart, the offset equal to the LTR permutation point. Merged AF 0.575, and the assemblies agree allele for allele.
- **chr6:78,894,316 polarity was inverted** — the masked reference holds a whole provirus.
- 80 non-HML-2 `LTR/ERVK` records were silently skipped and then passed `--strict` unfiltered; they are now labelled `other`.

### Validation

- Cluster runs on the CaG cohort at three revisions, all passing 41 assertions; the last exercises stage E end to end.
- Fixture selftests: 15 architecture cases, 16 consolidation assertions on real records.
- Consolidation validated against the graph calls at all five flagged loci; chr11 reproduces the assemblies exactly.
- Discovery outputs reproduce the existing published allele frequencies at all 22 unchanged loci.

### Caveats

- **Only the giraffe back end is validated** for consolidation. `graphaligner` and `pangenie` fail with an explicit error rather than an unchecked answer.
- **Tandem duplications carry no genotypes** by design (`--hervk_mask_tandem`), so those loci contribute `AN=0`.
- **`AN` can fall below `2N`** where a member was structurally uncalled. `AC` is unaffected — quote both, not `AF` alone.
- **`DENOVO_LTR` is reserved but not implemented.** Every degenerate CaG record is resolved by the reference check; it is a gap only for cohorts without a good reference annotation.
- **Sex-aware ploidy and PAR are untouched.** The reconciler parses whatever ploidy it is given and never imposes one; non-autosomal loci are flagged `PLOIDY_UNVERIFIED`.
- **`k` is an alignment property**, varies between haplotypes and callers, and nothing should key on its value.

### Please note

This branch also carries the **in-progress mkdocs restructure** (deletion of
`README_pipeline_full.md`, `docs/Makefile`, `docs/make.bat`, `docs/source/*`,
and renames under `docs/design-notes/`). Those were staged in the working tree
when the branch was cut and were swept into `9d19be0`. They are unrelated to
HERV-K and can be split out if you would rather they landed separately.
