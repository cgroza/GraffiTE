# HERV-K route coverage tests

`handout/` — the phase-1 discovery test. A GraffiTE `--human` run over an
existing PAV call set with `--genotype false`, checked against predictions
derived from the raw RepeatMasker tables. Self-contained; see
`handout/README_HPC_CLAUDE.md`.

Planned (not yet built):

- **Tier 1** — a ground-truthed synthetic mini-genome exercising svim-asm and
  PAV without reads.
- **Tier 2** — long reads: sniffles2 detection, GraphAligner and `lr-giraffe`
  genotyping. This is where the open sniffles questions get answered: whether
  its insertion ALTs are sequence-resolved enough for RepeatMasker to recover
  LTR termini at all, and how its breakpoint placement compares.
- **Tier 3** — short reads: PanGenie and `sr-giraffe`.

See the phase plan for the locus design and assertions.
