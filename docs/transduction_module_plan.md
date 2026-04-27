# Plan: TE Transduction Annotation Module for GraffiTE

## Background

**What is a 3' transduction?** When a source LINE-1 (or SVA, which uses L1 machinery) is
transcribed, the RNA polymerase sometimes reads past the element's own poly-A signal into the
downstream genomic flank before terminating. The resulting extended mRNA is reverse-transcribed by
TPRT, producing a new insertion that carries:

```
[TSD5'] [L1 / SVA body] [transduced genomic sequence] [poly-A tail] [TSD3']
```

The transduced sequence is genomically unique to the **source locus** — the specific reference L1
that spawned the retrotransposition event. This makes it both a biological signal (active source
elements, phylogenetic relationships among insertions) and an annotation problem amenable to
sequence mapping. Frequency: ~20% of full-length L1 retrotranspositions carry detectable 3'
transductions (Goodier et al. 2000; Tubio et al. 2014 *Science*; Gardner et al. 2017 *Science*).

**What about minus-strand insertions?** TPRT always uses the same L1 mRNA (5'→3' of the L1), but
the insertion can land on either strand of the target locus. For a **minus-strand** insertion the
ALT sequence (always written 5'→3' on the reference plus strand) is the reverse complement of the
mRNA order:

```
[polyT] [transduced_RC] [L1_RC]
```

The transduced sequence therefore appears at the **5' end** (low coordinates) of the ALT rather
than the 3' end.

**Existing signals in GraffiTE that enable this:** The pipeline already provides:
- `indels.fa` with full ALT sequences (one sequence per SV ID)
- Raw RepeatMasker `.out` with per-base masked coordinates *within* each ALT sequence
- `polyA=TRUE/FALSE` annotation distinguishing TPRT from other mechanisms <!--  yes but only happens once we have determined that combined RM+ULTRA coverage is >= threshold, which, in case of transduction may not be satisfied. We would need to evaluate poly-A on single-hit SV before filtering on RM+ULTRA coverage -->
- `L1_5PINV` for twin-priming / 5' inversions
- `matching_classes` identifying LINE/L1, Retroposon/SVA, SINE/Alu candidates
- The reference FASTA and minimap2 are already in the container

---

## Detection Algorithm

### Step 1 — Candidate selection

Select insertions satisfying all of:
- `SVTYPE=INS` <!-- It can also be that the reference genome has a transduction event that is absent from a sample we are comparing it to, and thus this transduction will appear as a deletion relative to the ref -->
- `matching_classes` contains `LINE/L1`, `Retroposon/SVA`, or `SINE/Alu`
- `polyA=TRUE` (confirms TPRT mechanism) <!--  yes but only happens once we have determined that combined RM+ULTRA coverage is >= threshold, which, in case of transduction may not be satisfied. We would need to evaluate poly-A on single-hit SV before filtering on RM+ULTRA coverage -->
- `abs(SVLEN) - total_match_length >= min_transduction_length` (enough unaccounted-for bases)

Alu transductions are possible but much rarer and shorter; treat them as candidates but expect low
yield.

### Step 2 — Parse RM coordinates within the ALT

From the raw RM `.out` file (collected across chromosomes), for each candidate SV ID:
1. Read all RM hits (columns: `qry_id`, `qry_start`, `qry_end`, `strand`)
2. Merge overlapping intervals → list of **masked intervals**
3. Compute **unmasked complement** against `[1, abs(SVLEN)]`
4. Determine the dominant TE strand (strand of the RM hit with the largest span)

### Step 3 — Locate the transduced sequence

| Insertion strand | Target interval |
|---|---|
| `+` | Last unmasked interval at the 3' end (highest coordinates), immediately before the poly-A tail |
| `−` | First unmasked interval at the 5' end (lowest coordinates), immediately after the polyT run |

Apply `min_transduction_length` cutoff (default 100 bp). If no qualifying unmasked interval
exists → not a candidate.

For **minus-strand** insertions, the polyT run is located by scanning low-coordinate unmasked
sequence for a T-run (>=10 nt) at the very start of the ALT; the transduced sequence is the
unmasked interval after it. <!--  we should also require that the transduced fragment is not a tandem repeat as identified by ULTRA; Simple_repeat and Low_complexity will show in the raw repeatmasker output though, but maybe we can use both RM and ULTRA as a safety -->

### Step 4 — Extract and map

1. Extract the candidate transduced subsequence from the ALT (VCF column 5)
2. Write to a temporary FASTA
3. Run `minimap2 --cs -c --secondary=no` against the reference
4. Parse the PAF output:
   - Filter: MAPQ >= `params.transduction_mapq` (default 20)
   - Filter: alignment identity >= 95% (NM / query_len from cs tag)
   - Exclude hits overlapping the insertion site ± 10 kb (self-mapping)
   - Keep the best hit (highest MAPQ, then longest alignment)
<!-- 
Here, ideally we want to do that from the original vcf before repeatmasker process as the transduced fragment can make repeat portion < threshold
   -->

### Step 5 — Optional source validation

If `params.ref_te_annotation` is provided (BED of reference TEs, e.g., from UCSC rmsk or a
dfam intersect):
- Intersect the source hit coordinates with reference TEs using bedtools
- Require overlap or proximity (<= `params.max_source_distance`, default 5 kb) to a reference
  L1 or SVA element
- Record the source element name/family in `TRANSD_SRC_TE`

### Step 6 — Annotate VCF and write report

Add INFO fields to the VCF via `bcftools annotate`; write a flat TSV transduction report.

---

## New INFO Fields

```
##INFO=<ID=TRANSDUCTION,Number=0,Type=Flag,Description="Insertion carries a 3-prime-transduced sequence detected by remapping to the reference">
##INFO=<ID=TRANSD_LEN,Number=1,Type=Integer,Description="Length (bp) of the transduced sequence">
##INFO=<ID=TRANSD_SRC_CHROM,Number=1,Type=String,Description="Chromosome of the inferred transduction source locus">
##INFO=<ID=TRANSD_SRC_START,Number=1,Type=Integer,Description="Start coordinate of the transduction source locus">
##INFO=<ID=TRANSD_SRC_END,Number=1,Type=Integer,Description="End coordinate of the transduction source locus">
##INFO=<ID=TRANSD_MAPQ,Number=1,Type=Integer,Description="Mapping quality of transduced sequence to source locus">
##INFO=<ID=TRANSD_SRC_TE,Number=1,Type=String,Description="TE annotation overlapping the transduction source locus (if --ref_te_annotation provided)">
```

---

## New Nextflow Process

```groovy
// in module/main.nf
process annotate_transductions {
  publishDir "${params.out}/3_TSD_search", mode: 'copy'

  input:
  path(pangenome_vcf)       // from concat_repeatmask.out.vcf_ch
  path("rm_out_*.out")      // collected raw RM .out files (new emit from repeatmask_VCF)
  path(ref_fasta)

  output:
  path("pangenome.transduction.vcf"), emit: vcf_ch
  path("transduction_report.tsv")

  script:
  def ref_te_arg = params.ref_te_annotation ? "--ref_te ${params.ref_te_annotation}" : ""
  """
  cat rm_out_*.out > all_rm.out
  detect_transductions.py \\
    --vcf ${pangenome_vcf} \\
    --rm_out all_rm.out \\
    --ref ${ref_fasta} \\
    --min_transd_len ${params.min_transduction_length} \\
    --mapq ${params.transduction_mapq} \\
    --max_source_dist ${params.max_source_distance} \\
    ${ref_te_arg} \\
    --out_vcf pangenome.transduction.vcf \\
    --out_tsv transduction_report.tsv
  """
}
```

### Required change to `repeatmask_VCF`

Add one new emit to expose the raw RM `.out` file (the filename is deterministic — RepeatMasker
always writes `<input>.out` inside the `-dir` directory):

```groovy
output:
  // ... existing outputs ...
  path("repeatmasker_dir/indels.fa.out"), emit: rm_out   // new
```

---

## New Script: `bin/detect_transductions.py`

Core logic outline (not a full implementation):

```python
# For each qualifying INS:
for sv_id, sv in candidates:
    hits     = rm_out[sv_id]                          # list of (start, end, strand)
    masked   = merge_intervals(hits)
    unmasked = complement(masked, 1, abs(svlen))
    strand   = dominant_strand(hits)                   # '+' or '-'

    if strand == '+':
        # transduced seq is at 3' end, before the poly-A tail
        interval = last_unmasked_before_polya(unmasked, polya_start)
    else:
        # transduced seq is at 5' end, after the polyT run
        interval = first_unmasked_after_polyt(unmasked, polyt_end)

    if interval is None or interval.length < min_transd_len:
        continue

    transd_seq = alt[interval.start : interval.end]
    paf_hits   = run_minimap2(transd_seq, ref_fasta)
    best_hit   = filter_and_rank(paf_hits, insertion_site, mapq_threshold)

    if best_hit:
        annotate_vcf(sv_id, TRANSDUCTION=True, TRANSD_LEN=interval.length, ...)
```

Key dependencies already in the GraffiTE container: `minimap2`, `bcftools`, `bedtools`, Python 3.

---

## Integration in `main.nf`

```groovy
// After concat_repeatmask, before genotyping:
if(params.transduction) {
  annotate_transductions(
    concat_repeatmask.out.vcf_ch,
    repeatmask_VCF.out.rm_out.collect(),
    ref_asm_ch
  )
  annotate_transductions.out.vcf_ch.set{ vcf_ch }
}
// vcf_ch then feeds into pangenie_index / make_graph as before
```

If `--graffite_vcf` is used (skipping discovery), the transduction module is not applicable
because the RM `.out` files are not produced. Document this limitation.

---

## New Parameters (`nextflow.config`)

```groovy
transduction             = false  // enable transduction annotation module
min_transduction_length  = 100    // min unmasked bp at 3'/5' end to attempt mapping
transduction_mapq        = 20     // min MAPQ for source locus hit
max_source_distance      = 5000   // max distance (bp) from source hit to a reference TE
ref_te_annotation        = false  // BED of reference TEs for source validation (optional)
```

---

## Output Files

| File | Location | Description |
|---|---|---|
| `pangenome.transduction.vcf` | `3_TSD_search/` | Full pangenome VCF with transduction INFO fields added |
| `transduction_report.tsv` | `3_TSD_search/` | Flat table: SV_ID, CHROM, POS, SVLEN, matching_classes, TRANSD_LEN, TRANSD_SRC_CHROM, TRANSD_SRC_START, TRANSD_SRC_END, TRANSD_MAPQ, TRANSD_SRC_TE |

---

## Caveats and Open Questions

### Biological

- **Truncated L1 insertions**: 5'-truncated L1s (the majority of non-reference insertions) may
  still carry a 3' transduction because the transduction comes from the source element, not the
  new copy. Detection should still work as long as `abs(SVLEN) - total_match_length` is large
  enough. However, short insertions (<500 bp) have less sequence context and higher false-positive
  risk.

- **SVA transductions**: SVA VNTRs are annotated as Simple_repeat by ULTRA/RM, so the unmasked
  tail after the VNTR is the correct target. The trusted-filter update (Simple_repeat bypasses
  ULTRA_TR_span threshold) is consistent with this — transduction detection should look at SVAs
  even with high ULTRA_TR_span.
  <!-- yes even though as I wrote before, we need to operate before filtering on % repeat threshold to be able to have candidates. Also, SVA-VNTR are unlikely to be loci with transduction: they are variation at VTNR only in fixed SVA so for SVA transduction, we actually want to look at other cases annotated as Retroposon/SVA-->

- **Alu transductions**: Alu uses L1 machinery and can transduce, but the events are shorter
  (Alu itself is only ~300 bp), rarer, and harder to distinguish from simple Alu dimers. Flag
  them but treat with lower confidence.

- **5' transductions**: ~2% frequency, more complex (often associated with twin-priming /
  L1_5PINV). Could be added later as a separate detection path; the `L1_5PINV` flag already
  provides a first-pass signal.


- **Somatic vs. germline**: GraffiTE is germline-focused. Cancer somatic transductions (Tubio
  2014) have very different allele frequencies and would need per-sample rather than
  population-level calling. Out of scope for now.

<!-- we'll dead with somatic later with a specific mode of GT-->

### Technical

- **RM `.out` filename stability**: `repmask_vcf.sh` always names the input `indels.fa` and
  passes `-dir repeatmasker_dir`, so RepeatMasker always writes `repeatmasker_dir/indels.fa.out`.
  The emit path is therefore deterministic.

- **minimap2 mode**: use `minimap2 -x asm5 --secondary=no --cs` for genomic sequence remapping.
  Avoid `map-ont` or `map-pb` — the transduced sequence is a short genomic fragment, not a long
  read.

- **Self-hit exclusion radius**: ±10 kb is conservative. If an active L1 is very close to its
  own insertion (uncommon but possible), this could exclude a valid hit. Could be made a
  parameter (`transduction_exclusion_radius`).

- **ID uniqueness**: relies on GraffiTE SV IDs being globally unique across chromosomes. They
  include genomic coordinates so they are, but `cat rm_out_*.out` joining per-chromosome files
  must not produce duplicate headers — the RM `.out` has 3-line headers per file; the
  `detect_transductions.py` parser should skip lines starting with `SW` or blank lines.

- **`--RM_dir` mode**: when users skip RepeatMasker with `--RM_dir`, the `.out` files must be
  present at `<RM_dir>/<chrom>/repeatmasker_dir/indels.fa.out`. The module should document this
  requirement and either collect from there or disable transduction annotation when `--RM_dir` is
  used without the `.out` files present.

### User added comments:

- I think the transduction detection should be handled in a separate, optional process, that interact with the repeatmasker one but do not respond to the filtering over total % SV coverage, otherwise we will miss the "cargo".
- You need to see what information can be created by unique process and what has to be separated, even if partially redundant. The fact that it'll be optional (or we can make it even also ONLY searching for transduction) will be known and accepted by the user. We can take advantage of nextflow to parallelize etc... 