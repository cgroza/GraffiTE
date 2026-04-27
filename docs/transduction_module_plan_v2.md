# Plan v2: TE Transduction Annotation Module for GraffiTE

Supersedes `transduction_module_plan.md` (v1) after a code re-review, integration of
user-supplied biological/technical comments, and a literature survey of prior
transduction-detection tools. Key reframing in v2:

- The transduction module **runs on a parallel branch** off `repeatmask_VCF`, on the
  **pre-coverage-filter** annotated VCF, so SVs whose transduced cargo drops
  `total_repeat_span` below the threshold are not lost.
- `add_polyA.py` already does the right thing on `n_hits==1` SVs — we just call it earlier
  in the transduction branch (it currently runs in `concat_repeatmask`, after coverage
  filtering).
- DEL SVs are treated symmetrically with INS (a transduction fixed in the reference can
  appear as DEL in another sample).
- The candidate transduced fragment must itself **not be a tandem repeat** (cross-checked
  against ULTRA + RM Simple_repeat/Low_complexity).
- The module is **fully optional** and **does not modify the main RM/TSD/concat path**;
  it only annotates the final pangenome VCF after the main path completes.
- A new **Section 12 ("Prior art")** maps the proposed algorithm against six published
  tools (TraFiC, TIGER, MELT-TRANSDUCTION, PALMER, TLDR, xTea), separates borrowed-and-
  credited recipe from GraffiTE-specific design choices, and lists the citations to use
  when describing the module in a future paper.

### Document layout

| § | Topic |
|---|---|
| 1 | What changed from v1 |
| 2 | Biology recap |
| 3 | Inputs already produced by `repeatmask_VCF` |
| 4 | Detection algorithm (revised, 6 steps) |
| 5 | New Nextflow processes |
| 6 | New Python script outline |
| 7 | `main.nf` integration |
| 8 | New parameters |
| 9 | New VCF INFO fields |
| 10 | Output files |
| 11 | Rejection-reason vocabulary |
| 12 | Prior art and citation guidance |
| 13 | Caveats and open questions |

---

## 1. What changed from v1

| v1 assumption | v2 fix | Reason |
|---|---|---|
| Run after `concat_repeatmask` on `pangenome.vcf` | Run **per-chromosome in parallel** with `tsd_search`, off `repeatmask_VCF` outputs | The pangenome VCF only contains SVs that passed `total_repeat_span > repeat_span_cutoff`; transduction cargo lowers that ratio so we need to evaluate candidates **before** the cutoff |
| `polyA=TRUE` is a usable filter | `polyA` must be **re-evaluated** by running `add_polyA.py` on the pre-filter VCF | `add_polyA.py` already runs in `concat_repeatmask` line 207, but only on the post-filter VCF; we run it earlier on the transduction branch, no algorithmic change |
| `SVTYPE=INS` only | `SVTYPE ∈ {INS, DEL}` | A reference-fixed transduction absent in another sample manifests as DEL relative to ref. `indels.fa` already contains both (REF for DEL, ALT for INS) |
| Locate unmasked tail and remap | Same, **plus** require the candidate fragment is not itself tandem-repeat | A VNTR/microsatellite tail will map to many genomic loci, generating false-positive sources |
| SVA candidates pass on `matching_classes` alone | SVA candidates pass only if the dominant RM hit is an SVA family (not Simple_repeat from a VNTR-only polymorphism) | SVA-VNTR polymorphisms are not transduction events; they are length variation in a fixed SVA's VNTR |
| `--graffite_vcf` simply unsupported | Same (documented), plus **`--RM_dir` is supported** because `genotypes_repmasked.vcf.gz` and `repeatmasker_dir/` are already published per chrom | `repeatmask_VCF` publishes everything to `2_Repeat_Filtering/<i>/`, so a re-run with `--RM_dir` has the inputs |
| One Python script handles everything | Same script, but called **per-chrom** to leverage Nextflow parallelism, then a small **merge** step | Aligns with the existing per-chrom RM scatter |

---

## 2. Biology recap (unchanged from v1, for context)

3' transduction during TPRT produces an insertion of the form:

```
  (+ strand insertion)   [TSD5'] [L1/SVA body 5'->3'] [transduced flank] [polyA] [TSD3']
  (- strand insertion)   [TSD5'] [polyT]    [transduced flank RC] [L1/SVA RC]   [TSD3']
                                  ^                                               ^
                                  alt[0..]                                        alt[end]
```

The transduced flank is genomically unique to the **donor locus** (the source L1 or SVA in
the reference). Detection = locate the unmasked tail in the ALT (or REF for DEL), remap to
the reference, call the source.

`add_polyA.py` already encodes the right window-anchored search: 8+ bp run within 5 bp of
the relevant terminus, ≥0.8 A/T purity, TSD-trimmed. Strand convention in
`RM_hit_strands`: `+` = polyA at 3' end, `C` = polyT at 5' end.

---

## 3. Inputs already available from `repeatmask_VCF`

The current `repeatmask_VCF` process (module/main.nf:232-258) already emits or publishes
everything the transduction path needs. **No new file generation required.**

| File | How produced | Currently exposed as |
|---|---|---|
| `genotypes_repmasked.vcf.gz` | `repmask_vcf.sh` lines 141-143 — full annotated VCF, **before** the `total_repeat_span` cutoff | `emit: repmasked_vcf_debug` (line 242) |
| `repeatmasker_dir/indels.fa` | `repmask_vcf.sh` lines 8-10 — both INS ALT and DEL REF sequences | inside `emit: vcf` (line 239) |
| `repeatmasker_dir/indels.fa.out` | RepeatMasker | inside `emit: vcf` |
| `ultra_out.bed` | `repmask_vcf.sh` lines 36-45 | `emit: ultra_bed` (line 240) |

**Recommended change**: bundle these four into one tuple emit named `transduction_in` so
the consumer process has a single, stable input shape:

```groovy
// in repeatmask_VCF (additive — does not break existing emits)
output:
  // ... existing emits unchanged ...
  tuple path("genotypes_repmasked.vcf.gz"),
        path("repeatmasker_dir/"),
        path("ultra_out.bed"),
        emit: transduction_in
```

The `repmasked_vcf_debug` emit can stay (or be aliased) for backward compatibility.

---

## 4. Detection algorithm (revised)

### Step 0 — Compute polyA on the pre-filter VCF

Run `add_polyA.py` on `genotypes_repmasked.vcf.gz` directly (no coverage filter). Output is
the same VCF with `polyA=TRUE/FALSE/NA` populated for every record where `n_hits==1`. This
is the only step that recomputes an existing annotation, and it's near-free.

### Step 1 — Candidate selection

Keep records satisfying **all** of:

```
SVTYPE in (INS, DEL)
n_hits == 1
matching_classes in (LINE/L1, Retroposon/SVA, SINE/Alu)
polyA == TRUE
abs(SVLEN) - total_match_length >= params.min_transduction_length
```

**SVA-specific guard**: for `matching_classes="Retroposon/SVA"`, also require that the
single RM hit's `repeat_id` is an SVA family (e.g. starts with `SVA`), not just a
Simple_repeat masquerading as one. In practice this is already enforced by `n_hits==1` +
the `matching_classes` filter (Simple_repeat hits are excluded from `matching_classes` in
`annotate_vcf.R`), but we assert it explicitly for safety.

**Alu**: kept as candidates with lower expected yield; flagged with `TRANSD_LOW_CONFIDENCE`
in the report.

### Step 2 — Parse RM .out and locate the candidate fragment

For each candidate SV ID:

1. Load all `repeatmasker_dir/indels.fa.out` rows for that ID (`qry_start`, `qry_end`,
   `strand`, `matching_class`, `repeat_id`).
2. Pick the **primary hit** = the row whose class matches the candidate's
   `matching_classes` and has the largest span. Read its strand.
3. Compute the **unmasked tail**:
   - `strand == '+'`: `tail = alt_seq[primary_hit.qry_end : len(alt_seq)]` (3' tail)
   - `strand == 'C'`: `head = alt_seq[0 : primary_hit.qry_start]` (5' head)
4. Locate the polyA/polyT inside the tail with the same logic as `add_polyA.py`
   (`MIN_LEN=8`, `MIN_PURITY=0.8`, `MAX_SLACK=5`). The transduced fragment is the slice
   between the primary hit and the polyA/polyT run:
   - `+`: `transd = alt_seq[primary_hit.qry_end : polyA.start]`
   - `C`: `transd = alt_seq[polyT.end : primary_hit.qry_start]`
5. Reject if `len(transd) < params.min_transduction_length`.

For DEL records the same logic applies — `indels.fa` keys DEL sequences (REF, line 10 of
`repmask_vcf.sh`) by SV ID exactly like INS. The "site" for self-hit exclusion (Step 4) is
the deletion span itself.

### Step 3 — Tandem-repeat exclusion of the candidate fragment

Reject the candidate if it is dominated by tandem repeat. Compute coverage of `transd` by:

- ULTRA intervals from `ultra_out.bed` for this SV ID, intersected with the fragment
  coordinates;
- RM Simple_repeat / Low_complexity hits from `indels.fa.out` (these classes are filtered
  out of `matching_classes` so they survive in the raw `.out` even though they don't
  contribute to `total_match_length`).

Reject if `(union TR coverage) / len(transd) > params.max_transd_tr_fraction` (default
0.5).

This is the v2 addition that addresses two of the user's comments at once: it prevents
SVA-VNTR-tail false positives, and prevents low-complexity tails from generating spurious
multi-mapping source calls.

### Step 4 — Remap to the reference

Write candidate fragments to a per-chrom FASTA:

```
>SV_ID|+|qstart-qend
ATCG...
```

Run minimap2 (already in container):

```
minimap2 -x asm5 -c --cs --secondary=no -t ${task.cpus} ref_fasta candidates.fa > candidates.paf
```

Filter the PAF:

- `MAPQ >= params.transduction_mapq` (default 20)
- alignment identity `>= params.transduction_min_identity` (default 0.95) computed from
  the `cs` tag (`matches / qlen`)
- exclude hits overlapping `[SV_POS - radius, SV_POS + |SVLEN| + radius]`, where
  `radius = params.transduction_exclusion_radius` (default 10000)
- keep the best remaining hit per candidate (highest MAPQ, tie-break on alignment length)

If no hit survives → record `TRANSDUCTION=NA` (candidate detected but unmapped). Optional:
emit a separate `transduction_orphans.tsv` for follow-up (could indicate a non-reference
source TE).

### Step 5 — Optional source TE validation

If `params.ref_te_annotation` (BED) is provided:

```
bedtools closest -d -t first -a best_hits.bed -b ${ref_te_annotation}
```

Record `TRANSD_SRC_TE=<repeat_name>` if distance ≤ `params.max_source_distance` (default
5000) and the matched element is in `{LINE/L1, Retroposon/SVA}`. Otherwise `TRANSD_SRC_TE=.`.

### Step 6 — Per-chrom output

Each `detect_transductions` task writes:

- `transduction_<chrom>.tsv` — one row per candidate (passed or rejected with reason)
- `transduction_annotations_<chrom>.tsv` — only the **passed** candidates, in the format
  consumed by `bcftools annotate` (CHROM, POS, ID, info kv pairs)

---

## 5. New processes

### `detect_transductions` (per chrom, optional, runs in parallel to `tsd_search`)

```groovy
process detect_transductions {
  publishDir "${params.out}/4_Transduction/per_chrom", mode: 'copy'

  input:
  tuple path("genotypes_repmasked.vcf.gz"),
        path("repeatmasker_dir/"),
        path("ultra_out.bed"),
        path(ref_fasta)

  output:
  path("transduction_*.tsv"),             emit: per_chrom_tsv
  path("transduction_annotations_*.tsv"), emit: per_chrom_annot

  script:
  def ref_te_arg = params.ref_te_annotation ? "--ref_te ${params.ref_te_annotation}" : ""
  """
  # 1) pre-filter polyA on the pre-coverage-filter VCF
  zcat genotypes_repmasked.vcf.gz > prefilter.vcf
  add_polyA.py prefilter.vcf -o prefilter_pa.vcf

  # 2) detect
  detect_transductions.py \\
    --vcf prefilter_pa.vcf \\
    --indels_fa repeatmasker_dir/indels.fa \\
    --rm_out repeatmasker_dir/indels.fa.out \\
    --ultra_bed ultra_out.bed \\
    --ref ${ref_fasta} \\
    --min_transd_len ${params.min_transduction_length} \\
    --max_tr_frac ${params.max_transd_tr_fraction} \\
    --mapq ${params.transduction_mapq} \\
    --min_identity ${params.transduction_min_identity} \\
    --exclusion_radius ${params.transduction_exclusion_radius} \\
    ${ref_te_arg} \\
    --threads ${task.cpus} \\
    --tsv_out transduction_${task.index}.tsv \\
    --annot_out transduction_annotations_${task.index}.tsv
  """
}
```

### `merge_transductions`

```groovy
process merge_transductions {
  publishDir "${params.out}/4_Transduction", mode: 'copy'

  input:
  path("transduction_*.tsv")
  path("transduction_annotations_*.tsv")

  output:
  path("transduction_report.tsv"),       emit: report
  path("transduction_annotations.tsv.gz"), emit: annot
  path("transduction_annotations.tsv.gz.tbi")

  script:
  """
  # full report
  awk 'NR==1 || FNR>1' transduction_*.tsv > transduction_report.tsv

  # annotation table: header + sorted body, bgzipped + tabixed for bcftools annotate
  awk 'NR==1 || FNR>1' transduction_annotations_*.tsv > transduction_annotations.tsv
  ( head -n1 transduction_annotations.tsv;
    tail -n +2 transduction_annotations.tsv | LC_ALL=C sort -k1,1 -k2,2n ) | \\
    bgzip > transduction_annotations.tsv.gz
  tabix -s1 -b2 -e2 -c'#' transduction_annotations.tsv.gz
  """
}
```

### `annotate_pangenome_transductions`

```groovy
process annotate_pangenome_transductions {
  publishDir "${params.out}/4_Transduction", mode: 'copy'

  input:
  path(pangenome_vcf)         // from concat_repeatmask.out.vcf_ch
  path(annot_gz)              // from merge_transductions.out.annot
  path(annot_tbi)

  output:
  path("pangenome.transduction.vcf"), emit: vcf_ch

  script:
  """
  # write header lines for the new INFO fields
  cat > transd_hdr.txt <<'EOF'
##INFO=<ID=TRANSDUCTION,Number=0,Type=Flag,Description="3' transduction detected by remapping the unmasked tail to the reference">
##INFO=<ID=TRANSD_LEN,Number=1,Type=Integer,Description="Length (bp) of the candidate transduced fragment">
##INFO=<ID=TRANSD_SRC_CHROM,Number=1,Type=String,Description="Chromosome of the inferred transduction source locus">
##INFO=<ID=TRANSD_SRC_START,Number=1,Type=Integer,Description="Start coordinate of the source locus (1-based)">
##INFO=<ID=TRANSD_SRC_END,Number=1,Type=Integer,Description="End coordinate of the source locus">
##INFO=<ID=TRANSD_MAPQ,Number=1,Type=Integer,Description="Mapping quality of the candidate to its source locus">
##INFO=<ID=TRANSD_IDENT,Number=1,Type=Float,Description="Alignment identity of the candidate to its source locus">
##INFO=<ID=TRANSD_SRC_TE,Number=1,Type=String,Description="Reference TE annotation overlapping or near the source locus, if --ref_te_annotation provided">
##INFO=<ID=TRANSD_LOW_CONF,Number=0,Type=Flag,Description="Low-confidence transduction call (e.g. Alu candidate, short fragment)">
EOF

  bcftools annotate \\
    -a ${annot_gz} \\
    -h transd_hdr.txt \\
    -c CHROM,POS,~ID,INFO/TRANSDUCTION,INFO/TRANSD_LEN,INFO/TRANSD_SRC_CHROM,INFO/TRANSD_SRC_START,INFO/TRANSD_SRC_END,INFO/TRANSD_MAPQ,INFO/TRANSD_IDENT,INFO/TRANSD_SRC_TE,INFO/TRANSD_LOW_CONF \\
    -Ov -o pangenome.transduction.vcf \\
    ${pangenome_vcf}
  """
}
```

---

## 6. New script: `bin/detect_transductions.py`

Outline (not full implementation):

```python
def main():
    args = parse_args()

    rm_hits   = parse_rm_out(args.rm_out)        # dict: sv_id -> [Hit(...)]
    ultra_iv  = parse_ultra_bed(args.ultra_bed)  # dict: sv_id -> [(s, e)]
    indels    = parse_fasta(args.indels_fa)      # dict: sv_id -> seq
    candidates = []

    for rec in read_vcf(args.vcf):
        if not passes_candidate_filter(rec, args):  # SVTYPE, n_hits, class, polyA, SVLEN
            continue
        sv_id = rec.id
        seq   = indels[sv_id]
        hit   = pick_primary_hit(rm_hits[sv_id], rec.matching_classes)
        frag  = locate_unmasked_tail(seq, hit, args)        # also confirms polyA in window
        if frag is None or len(frag.seq) < args.min_transd_len:
            log_reject(rec, 'no_tail_or_too_short'); continue
        if tr_fraction(frag, ultra_iv[sv_id], rm_hits[sv_id]) > args.max_tr_frac:
            log_reject(rec, 'tandem_repeat_fragment'); continue
        candidates.append((rec, frag))

    if not candidates:
        write_empty_outputs(args); return

    fasta = write_fasta(candidates, 'candidates.fa')
    paf   = run_minimap2(fasta, args.ref, args.threads)
    hits  = filter_paf(paf, candidates, args)               # MAPQ, identity, exclusion radius

    if args.ref_te:
        hits = annotate_source_te(hits, args.ref_te, args.max_source_dist)

    write_tsv(candidates, hits, args.tsv_out)
    write_annot_tsv(hits, args.annot_out)                   # CHROM POS ID + INFO kv
```

Existing helpers in repo to reuse: `add_polyA.has_anchored_tail` for polyA window detection.

---

## 7. main.nf integration

```groovy
// after concat_repeatmask
if(params.transduction) {
  // sanity check: transduction needs raw RM outputs from this run
  if(params.graffite_vcf) {
    error "--transduction is incompatible with --graffite_vcf (no raw RepeatMasker outputs)"
  }

  detect_transductions(repeatmask_VCF.out.transduction_in.combine(ref_asm_ch))

  merge_transductions(detect_transductions.out.per_chrom_tsv.collect(),
                      detect_transductions.out.per_chrom_annot.collect())

  annotate_pangenome_transductions(concat_repeatmask.out.vcf_ch,
                                   merge_transductions.out.annot,
                                   merge_transductions.out.annot.map{ it + '.tbi' })

  annotate_pangenome_transductions.out.vcf_ch.set{ vcf_ch }
}
```

For `--RM_dir` mode, `repeatmask_VCF` is skipped; we add a parallel branch that
reconstructs the transduction inputs from the published `${RM_dir}/<chrom>/`:

```groovy
if(params.RM_dir && params.transduction) {
  Channel.fromPath("${params.RM_dir}/*", type: "dir").map{ p ->
    tuple(file("${p}/genotypes_repmasked.vcf.gz", checkIfExists:true),
          file("${p}/repeatmasker_dir",            checkIfExists:true),
          file("${p}/ultra_out.bed",               checkIfExists:true))
  }.combine(ref_asm_ch).set{ transduction_in_ch }
  detect_transductions(transduction_in_ch)
  ...
}
```

(`genotypes_repmasked.vcf.gz` and `ultra_out.bed` are already published at
`${params.out}/2_Repeat_Filtering/<i>/` by the existing `publishDir` directive.)

---

## 8. New parameters (`nextflow.config`)

```groovy
// transduction module (off by default)
transduction                    = false   // enable transduction detection
transduction_only               = false   // optional: only run discovery + RM + transduction (skip TSD/concat/genotyping)
min_transduction_length         = 100     // min unmasked tail bp to attempt remapping
max_transd_tr_fraction          = 0.5     // max fraction of candidate fragment covered by ULTRA + RM Simple/Low_complexity
transduction_mapq               = 20      // min MAPQ for source locus hit
transduction_min_identity       = 0.95    // min alignment identity (matches / qlen)
transduction_exclusion_radius   = 10000   // bp around insertion site excluded from source search
ref_te_annotation               = false   // optional BED of reference TEs for source naming
max_source_distance             = 5000    // max bp from source hit to a reference TE
```

`transduction_only` short-circuits main.nf to skip `tsd_*`, `concat_repeatmask`, and
`pangenie/giraffe/graphaligner` blocks; only `repeatmask_VCF` + the three transduction
processes run, and `annotate_pangenome_transductions` writes its annotations onto a
**concatenated pre-filter VCF** instead of the trusted pangenome VCF.

---

## 9. New INFO fields on the final VCF

```
##INFO=<ID=TRANSDUCTION,Number=0,Type=Flag,...>
##INFO=<ID=TRANSD_LEN,Number=1,Type=Integer,...>
##INFO=<ID=TRANSD_SRC_CHROM,Number=1,Type=String,...>
##INFO=<ID=TRANSD_SRC_START,Number=1,Type=Integer,...>
##INFO=<ID=TRANSD_SRC_END,Number=1,Type=Integer,...>
##INFO=<ID=TRANSD_MAPQ,Number=1,Type=Integer,...>
##INFO=<ID=TRANSD_IDENT,Number=1,Type=Float,...>
##INFO=<ID=TRANSD_SRC_TE,Number=1,Type=String,...>
##INFO=<ID=TRANSD_LOW_CONF,Number=0,Type=Flag,...>
```

A candidate that passed all filters but failed to remap is **not** flagged
`TRANSDUCTION` in the VCF — it is recorded separately in `transduction_report.tsv`.

---

## 10. Output files

| File | Location | Description |
|---|---|---|
| `pangenome.transduction.vcf` | `4_Transduction/` | Final pangenome VCF with transduction INFO fields |
| `transduction_report.tsv` | `4_Transduction/` | All evaluated candidates with verdict (passed / rejected reason) |
| `transduction_annotations.tsv.gz(.tbi)` | `4_Transduction/` | bcftools-annotate-ready table of passed calls |
| `transduction_<i>.tsv` | `4_Transduction/per_chrom/` | Per-chrom raw output (for debugging) |

Downstream `vcf_ch` (consumed by `pangenie_index` / `make_graph` / `merge_VCFs`) is
overridden to the transduction-annotated VCF when `params.transduction` is set, so any
genotyping/graph step inherits the annotations transparently.

---

## 11. Rejection reasons logged in `transduction_report.tsv`

| Reason | Meaning |
|---|---|
| `not_candidate` | Failed Step 1 filter (multi-hit, wrong class, polyA=FALSE, SVLEN slack too small) |
| `no_primary_hit` | RM .out had no hit matching the declared class (should be rare) |
| `no_polyA_in_tail` | Strand-correct tail exists but no anchored A/T run (polyA recompute disagrees with annotation) |
| `tail_too_short` | Tail before the polyA < `min_transduction_length` |
| `tandem_repeat_fragment` | Tail dominated by ULTRA + RM low-complexity coverage |
| `unmapped` | minimap2 produced no hit |
| `low_quality_mapping` | All hits below MAPQ / identity threshold |
| `self_only` | Only hits were within the exclusion radius |
| `pass` | Candidate annotated as a transduction |

---

## 12. Prior art — published transduction-detection methods

Six pre-existing tools have already implemented L1/SVA transduction detection. The
GraffiTE module reuses a core algorithmic idea (remap the unmasked tail of an insertion to
the reference) that originated in analytical L1 surveys (Pickeral et al. 2000; Goodier et
al. 2000) and was first turned into NGS pipelines by TraFiC and TIGER. Credit should be
attributed when describing the method.

### Method comparison

| Tool | Reference | Input | Pop-scale germline | Somatic / cancer | Long-read | Short-read | 5' transduction | DEL-side detection | Source-element annotation | Tandem-repeat fragment guard |
|---|---|---|---|---|---|---|---|---|---|---|
| **TraFiC / TraFiC-mem** | Tubio et al. 2014 *Science* (+ Rodríguez-Martín et al. 2020 *Nat Genet* for PCAWG mem version) | tumor/normal BAM pairs | – | ✓ | – | ✓ | – | – | ✓ | – |
| **TIGER** | Tica et al. 2016 *BMC Genomics* | paired-end germline BAM | ✓ | – | – | ✓ | – | – | ✓ | – |
| **MELT-TRANSDUCTION** | Gardner et al. 2017 *Genome Res* | short-read BAM (1000G scale) | ✓ | – | – | ✓ | – | – | ✓ (FL-L1 source map) | – |
| **PALMER** | Zhou et al. 2020 *NAR* | long-read BAM (PacBio HiFi) | ✓ | – | ✓ | – | ✓ | – | ✓ (50 bp / 3.5 kb window) | – |
| **TLDR** | Ewing et al. 2020 *Mol Cell* | long-read BAM (ONT/PacBio) | ✓ | – | ✓ | – | ✓ | – | ✓ (`--trdcol` + `call_transductions.py`) | – |
| **xTea** | Chu et al. 2021 *Nat Commun* | short OR long read BAM | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ (realign clipped seq to FL-element flanks) | – |
| **GraffiTE-transduction (this plan)** | — | **pangenomic SV VCF** + RM `.out` + ULTRA `.bed` (assembled ALT/REF, no read access) | ✓ | (deferred) | (via assembled SVs) | (via assembled SVs) | (planned, not in v2) | **✓** | ✓ (optional `--ref_te_annotation`) | **✓** (ULTRA + RM Simple_repeat/Low_complexity) |

### What is borrowed vs. specific to GraffiTE

**Borrowed (and credited)** — the core "extract the unmasked tail bordered by polyA, remap
to the reference, call the source locus" recipe is shared with all six tools above. Within
that recipe:

- **Polyadenylation-anchored boundary detection** (find polyA at 3' / polyT at 5', the
  transduced flank is the slice between the TE and the tail) — universal in the field
  since Pickeral et al. 2000. `add_polyA.py` in GraffiTE already encodes this; we reuse it.
- **5'-transduction / minus-strand handling** (look for polyT at the 5' end of a
  reverse-complemented insertion) — explicit in PALMER, xTea, TLDR.
- **Source-element overlap with a reference TE BED** to assign donor identity — used by
  MELT, PALMER, xTea, TLDR.
- **Self-mapping exclusion radius** around the insertion site — used by TraFiC, TIGER,
  MELT.

**GraffiTE-specific design choices** (genuinely new in this module, not credited to
prior work):

- **VCF-first input model.** Every prior tool starts from BAMs and runs its own SV
  discovery. GraffiTE's module starts from an already-merged pangenomic SV VCF plus the
  pipeline's own RepeatMasker / ULTRA outputs. This is a structural shift, not an
  algorithmic improvement, but it is what lets GraffiTE annotate transductions without
  re-aligning reads.
- **Pre-coverage-filter evaluation.** Running `add_polyA.py` and the candidate selection
  on the un-cutoff VCF (`genotypes_repmasked.vcf.gz`) so that cargo-bearing insertions
  whose `total_repeat_span` falls below `repeat_span_cutoff` are still considered. Tools
  that operate on reads don't have this problem (they don't have an analogous global
  coverage filter).
- **DEL handling.** Treating a transduction fixed in the reference but absent from a
  sample as a deletion candidate. This is a natural consequence of the VCF-first input
  model — every prior tool was framed in terms of non-reference insertions.
- **ULTRA + RM Simple_repeat/Low_complexity union as a fragment-level tandem-repeat
  guard.** Prior tools rely on mapping uniqueness (MAPQ + best-hit) to suppress
  low-complexity false sources. GraffiTE adds an explicit upstream filter on the
  fragment itself, reusing the ULTRA bed and RM `.out` already produced by
  `repeatmask_VCF`. Inspired by xTea's "unique alignment" requirement but applied
  earlier in the pipeline.

### Suggested citations for the eventual GraffiTE paper

When describing the transduction module, cite:

1. **Biology / motivation**
   - Pickeral, Makalowski, Boguski, Boeke. 2000 *Genome Res* 10: 411-415 — first
     systematic description of L1-mediated 3' transduction.
   - Goodier, Ostertag, Kazazian. 2000 *Hum Mol Genet* 9: 653-657 — transduction
     frequency in L1 retrotransposition.
   - Macfarlane et al. 2013 *Hum Mutat* — transduction-driven SVA expansion.

2. **Algorithmic predecessors** (the prior tools whose recipe we re-implement)
   - Tubio et al. 2014 *Science* 345: 1251343 — TraFiC, first NGS-scale transduction
     pipeline (somatic).
   - Tica et al. 2016 *BMC Genomics* 17: 342 — TIGER, first germline computational
     transduction caller.
   - Gardner et al. 2017 *Genome Res* 27: 1916-1929 — MELT, with the
     MELT-TRANSDUCTION companion tool used in 1000 Genomes.
   - Zhou et al. 2020 *Nucleic Acids Res* 48: 1146-1163 — PALMER, long-read
     transduction detection on PacBio.
   - Ewing et al. 2020 *Mol Cell* 80: 915-928 — TLDR, long-read TE caller with
     transduction annotation.
   - Chu et al. 2021 *Nat Commun* 12: 3836 — xTea, cross-platform with explicit 3' and
     5' transduction modules.

3. **Reference catalogs (useful as comparison datasets, not as algorithmic sources)**
   - Rodríguez-Martín et al. 2020 *Nat Genet* 52: 306-319 — pan-cancer PCAWG
     transduction map (TraFiC-mem).
   - Niu et al. 2022 *Biology* 11: 1032 — 3' transduction map across 3202 human
     genomes.

A natural framing for the methods section: *"The GraffiTE transduction module
re-implements, in a pangenomic-VCF context, the polyA-anchored remapping strategy
established by TraFiC (Tubio et al. 2014) and TIGER (Tica et al. 2016) and extended to
long reads by PALMER (Zhou et al. 2020), TLDR (Ewing et al. 2020), and xTea (Chu et al.
2021). Distinct from these tools, our implementation operates downstream of SV calling on
an already-merged pangenomic VCF, evaluates candidates before the global repeat-coverage
filter so that cargo-bearing insertions are retained, and explicitly excludes
tandem-repeat fragments using the union of RepeatMasker low-complexity and ULTRA
annotations."*

---

## 13. Caveats and remaining open questions

### Biological

- **5' transductions** (~2% of L1 events, often associated with twin-priming /
  `L1_5PINV`) are not handled. A future extension: when `L1_5PINV` is set, also examine
  the unmasked region between the inverted segment and the 5' end.
- **Somatic transductions** (Tubio 2014 cancer L1s) need per-sample allele-frequency
  treatment. Out of scope (deferred to a future GraffiTE somatic mode).
- **Alu transductions**: kept as candidates with `TRANSD_LOW_CONF=1` because the events
  are short (~tens of bp) and easily confused with Alu dimers / read-through artifacts.
- **Truncated L1s**: 5'-truncated insertions can still carry a 3' transduction (the
  cargo comes from the source, not the new copy). The pre-coverage-filter design
  protects this — a 5'-truncated L1 with cargo can have `total_repeat_span` well below
  `repeat_span_cutoff` but is still a valid candidate.

### Technical

- **`indels.fa.out` parsing**: skip the 3-line RM header and any blank lines. Fields are
  whitespace-separated; the 1-based `qry_start`/`qry_end` columns must be converted to
  Python 0-based half-open before slicing the ALT string.
- **Strand convention mismatch**: RM .out uses `+`/`C`; ULTRA bed has no strand;
  `add_polyA.py` and `RM_hit_strands` both use `+`/`C`. The detection script should
  follow the same convention.
- **DEL site exclusion radius**: the "self" interval for a DEL is `[POS, POS+|SVLEN|]`
  rather than a point; expand the exclusion accordingly.
- **Per-chrom RM `.out` concatenation is no longer needed** in v2 (per-chrom processing
  by `detect_transductions`), so the v1 concern about duplicate headers from
  `cat rm_out_*.out` goes away.
- **minimap2 preset**: `-x asm5` (genomic, near-identical) is correct for typical
  transduction lengths (100 bp – ~5 kb). For very short fragments (<200 bp), MAPQ may
  be unreliable; the identity filter compensates.
- **`--graffite_vcf` incompatibility** is enforced with a hard error to avoid silent
  misuse.

### Future improvements (not in scope for first cut)

- Network-of-source-elements analysis: cluster candidates by source locus to identify
  the most active reference donors.
- 5' transduction detection (twin-priming aware).
- Per-sample transduction calls (somatic mode).
- Optional re-mapping with `splice` or `iter1` presets if `asm5` fails for unusual
  cases.
