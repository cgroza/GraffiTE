---
title: L1 5' inversions
description: Twin priming, the C/+ strand signature, and the L1_5PINV annotation.
---

# L1 5' inversions

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## Twin priming

Most L1 insertions are 5' truncated: reverse transcription starts at the polyA tail and stops
before reaching the 5' end of the RNA. In a fraction of them the 5' part of the element is also
inverted relative to the 3' part. Ostertag and Kazazian (2001) explained this by twin priming:
after the first nick, reverse transcription of the 3' end begins as usual, then the second
nicked strand primes a second reverse transcription on the same RNA, further upstream, in the
opposite direction. The two cDNAs are joined, and the element ends up as an inverted 5' piece
followed by a normally oriented 3' piece, often with a small deletion or duplication at the
junction between them.

<figure markdown="span">
![Twin priming: a second nick primes reverse transcription upstream on the same L1 RNA, in the opposite direction, so the integrated copy carries an inverted 5' segment ahead of the 3' segment](../assets/l1-twin-priming.png)
<figcaption>Twin priming. The 5′ segment is copied by a second priming event and integrates inverted.</figcaption>
</figure>

For a variant caller such an insertion is one event. For RepeatMasker it is two fragments of the
same L1 on opposite strands, and without a rule to join them the element reads as two hits.

## The strand signature

RepeatMasker reports each fragment's strand as `+` or `C`. An L1 with a 5' inversion appears as
a `C` fragment followed, along the variant sequence, by a `+` fragment, both belonging to the
same L1 (RepeatMasker gives them the same link ID). A `+` fragment followed by a `C` one is the
same element seen from the other strand: it is also a 5' inversion, but only the `C+` order is
recognised, because RepeatMasker reports the fragments in query order and the inverted piece
is 5' on the query.

<figure markdown="span">
![A RepeatMasker .out excerpt in which two fragments share one link ID, the first on strand C and the second on strand +](../assets/l1-strand-pattern.png)
<figcaption>The <code>C+</code> pattern in RepeatMasker output.</figcaption>
</figure>

## How GraffiTE detects it

`annotate_vcf.R` groups the fragments of each RepeatMasker link ID into one hit and records the
strands of its fragments, in query order, as one string (`C+`, `+`, `C`, `+C`). A hit is flagged
as an L1 5' inversion when its class is `LINE/L1` and that string is exactly `C+`. The rule does
not count the fragments: a `C` piece followed by two `+` pieces also qualifies.
<span class="src">`bin/annotate_vcf.R:80,95`</span>

The strand of the element as a whole is then inferred, because neither fragment's strand is the
answer. The consensus coordinate where each fragment begins is compared: if the `C` fragment
begins at a lower consensus position than the `+` fragment, the element is on the `+` strand;
otherwise it is on `C`. That inferred value is what `RM_hit_strands` reports for the hit, and
`polyA` is scanned on the end that strand implies. Every other hit reports its fragments'
strands as they came.
<span class="src">`bin/annotate_vcf.R:96-103`, `bin/add_polyA.py:107`</span>

<figure markdown="span">
![Consensus coordinates of the two fragments for a plus-strand and a minus-strand L1 with a 5' inversion](../assets/l1-consensus-coordinates.png)
<figcaption>Inferring the strand from where the two fragments sit on the L1 consensus.</figcaption>
</figure>

The rule runs on every dataset, not only human ones. It needs the library to name the class
`LINE/L1`, and it applies to any L1 subfamily in it.

## The L1_5PINV field

| Value | Meaning |
|---|---|
| `None` | No hit on the variant matched the rule (this is also the value for variants with no hit at all). |
| `<link ID>` | The RepeatMasker link ID of the hit flagged as inverted, the same number that appears in `RM_hit_IDs`, so the fragments can be found in `repeatmasker_dir/indels.fa.out`. Several IDs are comma-separated. |

<span class="src">`bin/annotate_vcf.R:147,171`, `bin/repmask_vcf.sh:134`</span>

Such a variant still has `n_hits=1` when the inverted L1 is its only element, so it is eligible
for the trusted and human subsets like any other single-hit L1. In v1.0 the same information was
written as `mam_filter_1=5P_INV`, and only with `--mammal`.
