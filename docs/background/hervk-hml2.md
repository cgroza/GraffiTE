---
title: HERV-K (HML-2) biology
description: Proviral and solo-LTR architecture, and how GraffiTE decides which allele a variant carries.
---

# HERV-K (HML-2) biology

!!! info "Applies to GraffiTE v1.1"
    HERV-K classification runs only under `--human`.

## HML-2 architecture

A HERV-K (HML-2) provirus is `LTR — INT — LTR`: two identical long terminal
repeats of 968 bp (`LTR5_Hs`; `LTR5A` is 1033 bp) flanking a 7,536 bp internal
region, 9,472 bp in all. Homologous recombination between the two LTRs excises
the intervening `LTR + INT` — exactly 8,504 bp, whatever the crossover point —
and leaves a single recombinant solo LTR.

So a locus sits in one of three states, and a structural variant call is always
a statement about **two** of them:

| state | length | |
|---|---|---|
| null | — | empty pre-integration site |
| solo LTR | 968 bp | one recombinant LTR |
| provirus | 9,472 bp | `LTR — INT — LTR` |

## Why size alone cannot tell them apart

An 8,504 bp insertion is equally consistent with `LTR+INT` entering a solo LTR
and with a complete provirus entering an empty site whose internal region
carries a deletion. GraffiTE therefore does not classify on size. It reads the
**architecture** of the variant allele out of the raw RepeatMasker fragments,
and where that is ambiguous it **masks the reference** at the locus and reads
off what is actually there.

Two architectures settle it outright:

- **Two full terminal LTRs** — both covering consensus 1–968. Nothing was
  consumed by the alignment, so the variant carries a whole provirus and the
  reference was empty. Reported as `ARCH_2LTR`.
- **One LTR split across the termini** — the 5′ fragment covering consensus
  `k+1..L` and the 3′ fragment `1..k`, complementary and summing to one LTR.
  The aligner broke a reference LTR at position `k`, so the reference held one.
  Reported as `ARCH_PERM`.

`HERVK_K` records `k`. It is a property of the **alignment**, not of the
biology: the inserted length is the same for every `k`, so the gap penalty
cancels and placement is decided by the handful of substitutions separating the
two LTR copies. It varies between haplotypes and between callers, and nothing
should key on its value.

At `k = 0` the split is degenerate — one whole LTR and nothing at the other end
— which is indistinguishable from a one-LTR-truncated provirus. About half of
real candidates land in that form, which is why the reference check is not a
fallback but the main mechanism.

## The classes

| `HERVK_CLASS` | REF ↔ ALT | |
|---|---|---|
| `null_solo` | null ↔ solo | a solo LTR segregates |
| `solo_prov` | solo ↔ provirus | dimorphic; no null observed |
| `truncated_prov` | as above, internal region incomplete | |
| `null_prov` | null ↔ provirus | a whole provirus segregates |
| `copy_number` | `prov_xN` ↔ `prov_xM` | the array gains or loses whole units |
| `other` | — | not HML-2, or unresolved |

`HERVK_ALLELE_REF` and `HERVK_ALLELE` give the states directly and are what
downstream analysis should read; `HERVK_CLASS` names the pair, not the
direction.

## Copy-number loci

An HML-2 locus is not always one provirus. chr7:4,699,540-4,717,514 (7p22.1a,
HERV-K108) carries two in tandem sharing a central LTR, and across 20 CaG
samples one, two and three units all segregate. `ref_n_units` counts the units
and `ref_unit_bp` measures the period, one internal region plus one LTR.

An SV whose `|SVLEN|` lands on a multiple of that period changes copy number by
that many units. Those states are written `prov_xN`, with `provirus` as the
N = 1 spelling. Zero units is a solo LTR, because an array of N units carries
N+1 LTRs and removing every unit leaves the one they shared, so the familiar
solo/provirus dimorphism is the N = 1 case of the same arithmetic.

The geometry, a unit of internal region plus exactly one LTR with the junction
inside the element the two copies share, is what unequal exchange between
misaligned units produces. Hughes and Coffin (2004, *PNAS* 101:1668-1672)
proposed that for the tandem HERV-K108 allele and suggested a solo LTR they
found was its reciprocal product. These calls belong in the analysis.

**Graph genotypes at these loci are withheld.** The ALT path repeats sequence
the reference already carries, so reads from the pre-existing copy traverse it
and non-carriers pick up ALT support. At chr6 the ALT fraction tracks provirus
dosage rather than carriage: provirus homozygotes run 0.13 to 0.42 while the
eight solo-LTR carriers sit at 0.00. Discovery genotypes come from
haplotype-resolved alignments, do not have this problem, and carry the allele
frequencies. Set `--hervk_mask_graph_gt_at_cnv false` to inspect the graph
calls anyway.

## Why the SVA hits are not SVA

RepeatMasker frequently reports an `SVA_*` hit at the terminus of an HML-2
element. This is homology, not a second element: SVA's SINE-R domain is
HERV-K LTR derived. Every such hit sits at SVA consensus ≥ ~900, and where it
is the sole annotation of a terminus its length completes the LTR exactly.
GraffiTE reassigns those to the LTR class and tiles hits winner-take-all on the
query axis, so an SVA hit lying on top of an `LTR5_Hs` hit is counted once.

## Library naming

The internal region is named `HERVK` by Dfam and `HERVK-int` by RepBase-derived
sets; `HERVK_int` and `HERVKint` also occur. All are recognised — **no modified
library is needed**. `HERVK9-int`, `HERVK11-int`, `HERVK14-int` and
`HERVKC4-int` are deliberately *not* matched: they are separate ERV lineages.
