---
title: Target site duplications
description: Why TSDs matter for mobile element insertions and how GraffiTE finds them.
---

# Target site duplications

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
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
<span class="src">`bin/TSD_Match_v2.sh:37-50`, `bin/exact_match.py:45-89,117`</span>

<figure>
--8<-- "assets/tsd-anatomy.svg"
<figcaption>
The two fragments, a duplication in the configuration that scores 0, and the score. In the
other snug configuration the two copies <em>start</em> at the junction, one at the beginning of
the variant and one at the beginning of the 3′ flank; it scores 1.
</figcaption>
</figure>

**`tsd_annotate_vcf.sh` writes `INFO/TSD`** for every record whose row in the summary says
`PASS`, as the 5' copy and the 3' copy comma-separated and upper-cased. Records that failed get
no `TSD` field.
<span class="src">`bin/tsd_annotate_vcf.sh:17-31`</span>

The search runs in batches of `--tsd_batch_size` variants per contig, and the per-batch
summaries and logs are concatenated into `3_TSD_search/TSD_summary.txt` and
`TSD_full_log.txt`.
<span class="src">`main.nf:157-162`, `module/main.nf:529-530`</span>

## Scoring and the PASS rule

Every match is scored by how far its two copies sit from the junction. With `R_start`, `R_end`,
`L_start`, `L_end` the 1-based columns of the match on each fragment:

```text
a     = ( |WIN - R_start| + |WIN - L_start| ) / 2
b     = ( |WIN - R_end|   + |WIN - L_end|   ) / 2
score = min(a, b)
```

`a` is low when both copies begin at the junction, `b` when both end there. The best candidate
for a variant is the lowest score, ties broken by the longest match. It passes when the score is
at most 5 bp. Both thresholds, and the 4 to 20 bp length range, are fixed in the script.
<span class="src">`bin/TSD_Match_v2.sh:48-50,86,112`</span>

A copy that ends exactly at the junction is at column `WIN`, so `b` is 0 for that configuration.
A copy that starts exactly at the junction is at column `WIN + 1`, so the start-anchored
configuration scores 1 rather than 0. The script's own comment records this asymmetry; the
threshold of 5 absorbs it.
<span class="src">`bin/TSD_Match_v2.sh:66-68`</span>

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
| 14 to 17 | `WIN` minus columns 8, 10, 9, 11: the start offsets of the R and L copies, then their end offsets |
| 18 | the score |
| 19, 20 | the sequence of the L copy and of the R copy |
| 21 | `PASS` or `FAIL` |

A variant with no exact match of 4 bp or more gets a shorter row: the ID, twelve `NA`, `no_hit`,
`no_hit`, `FAIL`. Read the file from the right (`$NF`, `$(NF-1)`, `$(NF-2)`) rather than by
column number, as `tsd_annotate_vcf.sh` does.
<span class="src">`bin/TSD_Match_v2.sh:50,57,112`, `bin/exact_match.py:92-109`</span>

`TSD_full_log.txt` shows, for each variant, both fragments over a base ruler, every candidate
match with its offsets, the chosen one, and the two fragments again with the copies underlined.
The column header printed there names 16 columns for rows that have 17, because it omits the bit
score; from the e-value on, read the header one column to the left.
<span class="src">`bin/TSD_Match_v2.sh:42,88-105`</span>

## Limitations

- **Exact matches only.** A duplication with a single mismatch between its copies, or one
  interrupted by a sequencing error in the assembly, is missed. v1.0 used an aligner here and
  reported partial matches; v1.1 traded that for a search that cannot produce spurious
  alignments in low-complexity flanks.
- **4 to 20 bp.** Shorter matches are below the seed; longer ones are treated as flank
  homology rather than a TSD.
- **One TSD per variant.** Only the best-scoring candidate is reported, so a variant with two
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
  <span class="src">`bin/add_polyA.py:74-85`</span>
