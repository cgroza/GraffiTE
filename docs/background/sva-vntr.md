---
title: SVA VNTR polymorphisms
description: Why VNTR-only SVA variants are reclassified, and the per-subfamily VNTR coordinates.
---

# SVA VNTR polymorphisms

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `ee7da10`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## SVA structure

An SVA element is a composite: a hexamer repeat at the 5' end, an *Alu*-like region, a
variable number of tandem repeats (the VNTR), a SINE-R region derived from an endogenous
retrovirus, and a polyA tail. The VNTR is made of GC-rich units of roughly 37 to 50 bp, and its
copy number varies between elements and between alleles of the same element.

<figure markdown="span">
![Structure of an SVA element: hexamer repeat, Alu-like region, VNTR, SINE-R, polyA tail](../assets/sva-vntr-structure.png)
<figcaption>The regions of an SVA element. The VNTR is the part that expands and contracts.</figcaption>
</figure>

## VNTR-only polymorphisms

Because the VNTR is a tandem array, it changes length by replication slippage or unequal
recombination without any transposition. Between two haplotypes that both carry the same SVA,
one can have a longer VNTR than the other, and an SV caller reports that difference as an
insertion or a deletion of a few hundred bp. RepeatMasker then annotates the variant sequence
as SVA, since VNTR sequence is SVA sequence, and without a rule to catch it the variant counts as
an SVA insertion polymorphism when no element moved.

GraffiTE recognises these by where the hit lands on the consensus: a variant whose single SVA
fragment sits entirely inside the VNTR region of its subfamily is a VNTR-only polymorphism.

## VNTR coordinates by subfamily

The VNTR interval on each subfamily consensus, as fixed in `annotate_vcf.R`:

| Subfamily | VNTR start | VNTR end |
|---|---|---|
| `SVA_A` | 436 | 855 |
| `SVA_B` | 431 | 867 |
| `SVA_C` | 432 | 851 |
| `SVA_D` | 432 | 689 |
| `SVA_E` | 428 | 864 |
| `SVA_F` | 435 | 857 |

<span class="src">`bin/annotate_vcf.R:132-134`</span>

The coordinates are those of the Dfam human consensus sequences with those names, so the rule
only fires for a library that uses them. A hit named otherwise (`SVA_F1`, or a species-specific
name) is left as is.

## How GraffiTE reclassifies them

A hit is a VNTR-only polymorphism when all of these hold:
<span class="src">`bin/annotate_vcf.R:136-149`</span>

- its class is `Retroposon/SVA`;
- it is a single RepeatMasker fragment;
- its name is one of the six subfamilies above;
- its start on the consensus is greater than the subfamily's VNTR start, and its end is less
  than the VNTR end (strictly inside on both sides).

Such a hit is then written with a `(VNTR_only)` suffix on its name in `repeat_ids`, and its
entry in `matching_classes` becomes `Simple_repeat` instead of `Retroposon/SVA`. The hit is
still counted in `n_hits` and still contributes to `total_repeat_span`.
<span class="src">`bin/annotate_vcf.R:151-156`</span>

```text
repeat_ids=SVA_F(VNTR_only);matching_classes=Simple_repeat
```

That relabelling is what the downstream filters see:

- The **trusted subset** admits a `Simple_repeat` record without testing its `ULTRA_TR_span`,
  so a VNTR expansion (which is all tandem repeat) is not excluded on that ground.
  <span class="src">`module/main.nf:14`</span>
- The **`--human` filter** admits `Simple_repeat` records whose `repeat_ids` match
  `--human_sva_ids` (default `^SVA_[DEF]`), without requiring a polyA tail, so VNTR
  polymorphisms of the young subfamilies reach `pangenome.human.vcf`, distinguishable from
  SVA insertions by the suffix.
  <span class="src">`module/main.nf:520,522,528`, `nextflow.config:67`</span>

To count SVA insertion polymorphisms alone, exclude them:

```bash
bcftools view -e 'INFO/repeat_ids~"VNTR_only"' pangenome.human.vcf
```

## The subsets truncate this set

!!! warning "`pangenome.trusted.vcf` and `pangenome.human.vcf` are not a VNTR catalogue"
    On the dataset below they keep one VNTR-only record in eleven, and the ones they keep are
    the longest. Read `pangenome.vcf` for VNTR work. A later release will size these records in
    VNTR units instead of bp.

Both subsets require `|SVLEN|` of at least 250 bp (`--trusted_min_svlen`,
`--human_min_svlen`). We set that threshold for insertions, where 250 bp is a truncated *Alu*.
The VNTR unit is about 49 bp, so the same threshold asks a VNTR change to span five units, and
most changes are smaller.

Measured on the 20-genome HPRC assembly set (PAV discovery, graph genotyping) used for the v1.1
human analysis:

| | records |
|---|---|
| `repeat_ids~"VNTR_only"` in the annotated callset | 1,141 |
| of those, `\|SVLEN\| < 250` bp, so dropped by the default | **1,039 (91%)** |
| median length change | 126 bp |

The length distribution peaks at 50, 78, 98, 127, 176 and 225 bp. Those are sums of a 49 bp and
a 78 bp VNTR unit, with 127 = 49 + 78 as the higher-order repeat, and the 250 bp default sits
above every one of them.

`--human_sva_ids` cuts further. GraffiTE matches it against `repeat_ids`, which for a VNTR-only
record names the consensus the VNTR sequence scored best against. That is not the subfamily of
the element hosting it: on the dataset above the two agree for 402 of 1,138 records (35%), and
65% of the VNTR changes sitting in an SVA_A, SVA_B or SVA_C element are annotated to an old
consensus. The default `^SVA_[DEF]` therefore drops them. That filter alone takes SVA_A host
elements from 20.4% carrying a variable VNTR to 3.2%, which looks like a difference between
subfamilies and is not one.

With both filters off, 20% of reference SVA_A elements of at least 1 kb in CHM13v2 carry a
variable VNTR across the 20 samples, and 54% of SVA_F elements do. With both filters on, 0% and
10%.

A floor remains below that, from the SV caller. The PAV path keeps `|SVLEN| > 50`
<span class="src">`module/main.nf:161`</span>, and the sniffles and svim-asm paths ask for 100 bp
<span class="src">`module/main.nf:105,122,178`</span>. On the run above the smallest annotated
VNTR record is 50 bp, about one unit, so single-unit changes sit at the edge of what the
callset records.

### Working with the full set

Take the VNTR-only records from `pangenome.vcf`, which carries no size or subfamily filter, and
apply your own rules:

```bash
bcftools view -i 'INFO/repeat_ids~"VNTR_only"' pangenome.vcf
```

Genotypes for the same records are in the merged genotyped VCF, which GraffiTE annotates from
`pangenome.vcf` and likewise does not filter.

To get them inside `pangenome.human.vcf` instead, lower the threshold and widen the subfamily
list:

```bash
nextflow run cgroza/GraffiTE --human \
  --human_min_svlen 50 \
  --human_sva_ids '^SVA_'
```

Both changes loosen the filter for SVA insertions too, which is usually not what you want, so
prefer reading `pangenome.vcf` unless you need the whole subset relaxed.

In v1.0 the same records were flagged with `mam_filter_2=VNTR_ONLY:<family>:<start>:<end>`,
only with `--mammal`, and kept their `Retroposon/SVA` class.
