---
title: Target site duplications
description: Why TSDs matter for mobile element insertions and how GraffiTE finds them.
---

# Target site duplications

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `cc1f3ac`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## What a TSD is

When a non-LTR retrotransposon integrates by target-primed reverse transcription, the
endonuclease nicks the two strands of the target site a few bases apart. Repairing that offset
copies the bases between the nicks, so the new element ends up flanked by two identical short
direct repeats, the target site duplication. LTR retrotransposons and DNA transposons make one
too, through their own integrases, with a length characteristic of the family. A TSD is
therefore evidence that a sequence got where it is by transposition, rather than by any other
kind of structural change.

In a variant call, the two copies straddle the breakpoint. The reference genome holds one copy
(the original target site), and the inserted sequence carries the other at one of its ends,
depending on where the caller placed the breakpoint. GraffiTE looks for that pair.

## The search procedure

The search runs on every record that passed the repeat-span filter, whatever its `n_hits`, and
uses only exact matches. Three scripts do the work, one process each.

**`prepTSD.sh` builds two FASTA files.** For each variant, `tsd_flanks.py` reads `--tsd_win` bp
of reference ending at `POS` (the anchor base included) and `--tsd_win` bp starting just after
the variant: after `POS` for an insertion, after the deleted interval for a deletion. A window
that runs off the start of a contig is clamped and comes out shorter. Separately, the variant
sequence itself (ALT for an insertion, REF for a deletion) is trimmed to its first and last
`--tsd_win` bp; a variant no longer than the window is used whole. The default window is 30 bp.
<span class="src">`bin/prepTSD.sh:43-56`, `bin/tsd_flanks.py:50-59`, `nextflow.config:49`</span>

**`TSD_Match_v2.sh` compares two fragments per variant.** The L fragment is the 5' flank followed
by the first bases of the variant; the R fragment is the last bases of the variant followed by
the 3' flank. Each is two windows long and the junction between flank and variant sits at column
`WIN` in both. `exact_match.py` then lists every maximal exact match of at least 4 bp between R
and L on the plus strand, in BLAST tabular form, and matches longer than 20 bp are discarded.
<span class="src">`bin/TSD_Match_v2.sh:37-69`, `bin/exact_match.py:45-89,117`</span>

<figure>
--8<-- "assets/tsd-anatomy.svg"
<figcaption>
The two fragments, a duplication in the configuration that scores 0, and the score. In the
other snug configuration the two copies <em>start</em> at the junction, one at the beginning of
the variant and one at the beginning of the 3′ flank; it scores 0 too.
</figcaption>
</figure>

**`tsd_annotate_vcf.sh` writes `INFO/TSD`** for every record whose row in the summary says
`PASS`, as the 5' copy and the 3' copy comma-separated and upper-cased. Records that failed get
no `TSD` field. The header declares `Number=2`, one value per copy. `fix_vcf.py` rewrites
`pangenome.vcf` afterwards with vcfpy, which percent-encodes a comma inside a `Number=1` value.
In output from before commit `d1dfd7b` the header said `Number=1`, and the field shows as
`GATTACAG%2CGATTACAG`.
<span class="src">`bin/tsd_annotate_vcf.sh:22-36`</span>

The search runs in batches of `--tsd_batch_size` variants per contig, and the per-batch
summaries and logs are concatenated into `3_TSD_search/TSD_summary.txt` and
`TSD_full_log.txt`.
<span class="src">`main.nf:149-154`, `module/main.nf:573-574`</span>

## Scoring and the PASS rule

Every match is scored by how far its two copies sit from the junction. With `R_start`, `R_end`,
`L_start`, `L_end` the 1-based columns of the match on each fragment:

```text
a     = ( |WIN + 1 - R_start| + |WIN + 1 - L_start| ) / 2
b     = ( |WIN - R_end|       + |WIN - L_end|       ) / 2
score = min(a, b)
```

`a` is low when both copies begin at the junction, `b` when both end there. A candidate shorter
than 6 bp competes only if it scores 0.5 or less. Among the rest, candidates scoring 1.5 or less
count as ties, and the longest of them is the best. When none scores that low, the best is the
lowest score, ties broken by the longest match. The best candidate passes when its score is at
most 5 bp. These thresholds, and the 4 to 20 bp length range, are fixed in the script.
<span class="src">`bin/TSD_Match_v2.sh:48-49,59-60,69,115-116,142`</span>

A copy that ends exactly at the junction is at column `WIN`, and one that starts exactly there is
at column `WIN + 1`, so both snug configurations score 0. The tie at 1.5 is for insertions that
end in a poly(A) tail. When the TSD also opens with As, the tail and the TSD run together and the
boundary between them is uncertain by a base or two. A 4 bp run of A scoring 0.5 lower would
otherwise beat the real copy. At 1.5 the two copies can sit up to three bases off the junction
between them.
<span class="src">`bin/TSD_Match_v2.sh:61-68,104-116`, `test/tsd/test_tsd_match.sh`</span>

The rule for short candidates comes from a null test on the CaG set (20 genomes, CHM13v2.0). We
replaced each insertion's R window with that of another insertion of the same family, strand and
SV type, so that no duplication could span the pair, and ran the search unchanged. Without the
rule it passed a match in 72% of swapped Alu, L1 and SVA pairs, most of them 4 or 5 bp long and
off the junction, while real calls of 9 bp or more scored 0 or 0.5 in 97% of cases. With it,
36% of swapped pairs pass, and 53 calls on 5,614 real Alu, L1 and SVA records are lost, all
short and off the junction. Swapped pairs still get a call of 9 bp or more 1.6% of the time, so a
long TSD is good evidence of transposition and a short one is weak even when it sits on the
junction.
<span class="src">`bin/TSD_Match_v2.sh:50-60,115-116,142`, `test/tsd/test_tsd_match.sh`</span>

## Reading TSD_summary.txt

One row per variant searched. A row with a hit has 21 tab-separated columns:

| Column | Content |
|---|---|
| 1 | variant ID |
| 2, 3 | `R\|3P_end`, `L\|5P_end` (query and target fragment names) |
| 4 | identity, always `100.000` |
| 5 | match length in bp |
| 6, 7 | mismatches and gap opens, always `0` |
| 8, 9 | start and end of the copy on the R fragment |
| 10, 11 | start and end of the copy on the L fragment |
| 12, 13 | placeholder e-value and bit score, derived from the length |
| 14 to 17 | the start offsets of the R and L copies (`WIN + 1` minus columns 8 and 10), then their end offsets (`WIN` minus columns 9 and 11) |
| 18 | the score |
| 19, 20 | the sequence of the L copy and of the R copy |
| 21 | `PASS` or `FAIL` |

A variant with no exact match of 4 bp or more gets a shorter row: the ID, twelve `NA`, `no_hit`,
`no_hit`, `FAIL`. Read the file from the right (`$NF`, `$(NF-1)`, `$(NF-2)`) rather than by
column number, as `tsd_annotate_vcf.sh` does.
<span class="src">`bin/TSD_Match_v2.sh:69,76,142`, `bin/exact_match.py:92-109`</span>

`TSD_full_log.txt` shows, for each variant, both fragments over a base ruler, every candidate
match with its offsets, the chosen one, and the two fragments again with the copies underlined.
The column header printed there names 16 columns for rows that have 17, because it omits the bit
score. The header lines up through the e-value in column 11; column 12 of a row is the unnamed
bit score; from column 13 on, the right name is one header column to the left.
<span class="src">`bin/TSD_Match_v2.sh:42,118-135`</span>

## Limitations

- **Exact matches only.** A duplication with a single mismatch between its copies, or one
  interrupted by a sequencing error in the assembly, is missed. v1.0 used an aligner here and
  reported partial matches; v1.1 traded that for a search that cannot produce spurious
  alignments in low-complexity flanks.
- **4 to 20 bp.** Shorter matches are below the seed; longer ones are treated as flank
  homology rather than a TSD. A match under 6 bp counts only when it sits on the junction.
- **One TSD per variant.** Only the best candidate is reported, so a variant with two
  plausible duplications shows one.
- **The window bounds what can be seen.** A copy further than `--tsd_win` from the breakpoint,
  or a breakpoint the caller placed more than 5 bp from the true junction, does not pass. Raising
  `--tsd_win` widens the search and, since the score is measured from the junction, does not
  change what a snug duplication scores.
- **The variant sequence must be resolved.** A symbolic `<INS>` has no ends to compare;
  Stage A drops those.
- **PolyA tails can hide a copy.** The 3' end of a plus-strand Alu or L1 is a run of A; a TSD
  rich in A can be found at the wrong offset inside it. `polyA` is computed after the `TSD` copy
  is trimmed, in the other direction, for this reason.
  <span class="src">`bin/add_polyA.py:75-86`</span>
