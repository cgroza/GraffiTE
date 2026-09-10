---
title: HERV-K (HML-2) biology
description: >-
  Proviral and solo-LTR architecture, copy-number arrays, and how GraffiTE decides which two
  states of a locus a variant call compares.
---

# HERV-K (HML-2) biology

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).
    HERV-K classification runs only under `--human`; the procedure is in
    [Human MEIs](../guides/human-mei.md).

## HML-2 architecture

A HERV-K (HML-2) provirus is `LTR - INT - LTR`: two identical long terminal repeats of 968 bp
(`LTR5_Hs`; `LTR5A` is 1,033 bp) flanking a 7,536 bp internal region, 9,472 bp in all
<span class="src">`bin/hervk_arch.py:56-57`, `bin/hervk_classify.py:49-51`</span>. Homologous
recombination between the two LTRs excises the intervening `LTR + INT`, 8,504 bp whatever the
crossover point, and leaves a single recombinant solo LTR.

<figure>
--8<-- "assets/hervk-locus-states.svg"
<figcaption>
The states one locus can hold, and the value each carries in <code>HERVK_ALLELE_REF</code> and
<code>HERVK_ALLELE</code>. A fifth value, <code>partial</code>, names a reference element that is
neither a whole LTR nor a provirus with both LTRs.
</figcaption>
</figure>

A structural variant call is always a statement about **two** of these states, and the
polymorphism type depends on which two. A `null` allele against anything else is an insertion
polymorphism, the event behind the difference being a transposition; that is what the
`HERVK_MEI` flag records, and what makes a HERV-K locus comparable to an Alu or L1 insertion.
Solo LTR against provirus, or one unit against two, is structural variation in an element every
haplotype carries <span class="src">`bin/hervk_reconcile.py:136-170`</span>.

## Why size alone cannot tell them apart

An 8,504 bp insertion is equally consistent with `LTR + INT` entering a solo LTR and with a
complete provirus entering an empty site whose internal region carries a deletion. GraffiTE
therefore does not classify on size. It reads the **architecture** of the variant allele out of
the raw RepeatMasker fragments, and where that is ambiguous it **masks the reference** at the
locus and reads off what is there.

Two architectures settle it outright <span class="src">`bin/hervk_arch.py:254-356`</span>:

- **Two full terminal LTRs**, both covering consensus 1 to 968. Nothing was consumed by the
  alignment, so the variant carries a whole provirus. Reported as `ARCH_2LTR`.
- **One LTR split across the termini**: the 5′ fragment covering consensus `k+1..L` and the 3′
  fragment `1..k`, complementary and summing to one LTR. The aligner broke a reference LTR at
  position `k`, so the reference held one. Reported as `ARCH_PERM`. The same split inside the
  internal region is `ARCH_INT_PERM`, with its point in `HERVK_J`.

Which side of the pair the variant carries depends on its polarity. For an insertion the
architecture describes the allele being added; for a deletion it describes reference sequence
being removed, so `ARCH_2LTR` reads `null -> provirus` on an `INS` and `provirus -> null` on a
`DEL` <span class="src">`bin/hervk_classify.py:322-363`</span>.

`HERVK_K` records `k`. It is a property of the **alignment**, not of the biology: the inserted
length is the same for every `k`, so the gap penalty cancels and the handful of substitutions
separating the two LTR copies decides the placement. It varies between haplotypes and
between callers, and nothing should key on its value.

At `k = 0` the split is degenerate: one whole LTR and nothing at the other end, which is
indistinguishable from a one-LTR-truncated provirus. The reference check is therefore the main
mechanism, not a fallback: the classifier masks the reference in a window of `--hervk_ref_flank`
(1,500 bp) either side of the candidate, widens it when an element runs into the edge, and calls
the state from the LTR and internal-region base pairs the window holds
<span class="src">`bin/hervk_ref_state.py:46-72, 272-292`</span>.

## The classes

`HERVK_ALLELE_REF` and `HERVK_ALLELE` give the states directly and are what downstream analysis
should read. `HERVK_CLASS` names the pair, not the direction
<span class="src">`bin/hervk_classify.py:274-288, 424-467`</span>:

| `HERVK_CLASS` | states | |
|---|---|---|
| `null_solo` | null, solo | a solo LTR segregates |
| `solo_prov` | solo, provirus | dimorphic; no null observed |
| `truncated_prov` | as `null_prov` or `solo_prov`, with less than 80 % of the internal region | |
| `null_prov` | null, provirus | a whole provirus segregates |
| `copy_number` | `prov_xN`, `prov_xM` | the array gains or loses whole units |
| `other` | | not resolved to two HML-2 states |

`HERVK_EVIDENCE` says what settled the call: `ARCH_2LTR`, `ARCH_PERM`, `ARCH_INT_PERM`,
`ARCH_SOLO`, `CNV_PERIOD`, `REF_ANNOT` (the reference decided), `UNRESOLVED`, or `NON_HML2` (a
candidate with fewer than 50 bp of HML-2 sequence) <span class="src">`bin/hervk_classify.py:434-436`</span>.
`HERVK_PMAP` is a confidence under a size model and never decides a class.

## Copy-number loci

An HML-2 locus is not always one provirus. In the CaG cohort, the locus at chr7:4.70 Mb carries
two units in tandem sharing a central LTR, and one, two and three units all segregate
<span class="src">`module/main.nf:509-516`</span>. The reference step counts the units of a
proviral reference element and measures their period, one internal region plus one LTR
(`HERVK_N_UNITS_REF`, `HERVK_UNIT_BP`) <span class="src">`bin/hervk_ref_state.py:238-270`</span>.

An SV whose `|SVLEN|` lands on a whole number of that period, within 50 bp or 0.5 % of its
length, changes copy number by that many units. Those states are written `prov_xN`, with
`provirus` as the `N = 1` spelling. Zero units is a solo LTR, because an array of N units carries
N + 1 LTRs and removing every unit leaves the one they shared, so the familiar solo/provirus
dimorphism is the `N = 1` case of the same arithmetic
<span class="src">`bin/hervk_classify.py:225-271`</span>.

The geometry, one internal region plus exactly one LTR with the junction inside the LTR the two
copies share, is what unequal exchange between misaligned units produces. Hughes and Coffin
(2004, *PNAS* 101:1668-1672) proposed that mechanism for the tandem HERV-K108 allele and
suggested a solo LTR they found was its reciprocal product.

**Graph genotypes at these loci are withheld** (`HERVK_GT_MASKED`). The ALT path repeats
sequence the reference already carries, so reads from the pre-existing copy traverse it and
non-carriers pick up ALT support; in the CaG run the ALT fraction at such a locus tracked
provirus dosage rather than carriage. Discovery genotypes come from haplotype-resolved
alignments, do not have this problem, and carry the allele frequencies
(`HERVK_AC_DISC`, `HERVK_AN_DISC`). Set `--hervk_mask_graph_gt_at_cnv false` to keep the graph
calls <span class="src">`nextflow.config:98-106`</span>.

## Why the SVA hits are not SVA

RepeatMasker frequently reports an `SVA_*` hit at the terminus of an HML-2 element. This is
homology, not a second element: SVA's SINE-R domain is HERV-K LTR derived. GraffiTE treats an SVA hit
that starts at SVA consensus position 900 or later and abuts HML-2 sequence within 50 bp as
SINE-R, reassigns it to the LTR class and tiles hits
winner-take-all along the variant, so an SVA hit lying on top of an `LTR5_Hs` hit is counted
once <span class="src">`bin/hervk_arch.py:69-84, 201-233`</span>. The same homology is why the
`--human` filter carries its HERVK + SVA carve-out ([Human MEIs](../guides/human-mei.md)).

## Library naming

The internal region is named `HERVK` by Dfam and `HERVK-int` by RepBase-derived sets;
`HERVK_int` and `HERVKint` also occur. The classifier recognises all of them, so no modified
library is needed <span class="src">`bin/hervk_arch.py:65`</span>. `HERVK9-int`, `HERVK11-int`,
`HERVK14-int` and `HERVKC4-int` are deliberately *not* matched: they are separate ERV lineages.
The LTR families recognised are `LTR5_Hs`, `LTR5A`, `LTR5B` and `LTR5`
<span class="src">`bin/hervk_arch.py:56`</span>.

One place in the pipeline is stricter than the classifier: the carve-out in the `--human` filter
requires the literal name `HERVK-int`; see the note in [Human MEIs](../guides/human-mei.md).
