# GraffiTE documentation restructure — proposal

**Status:** Draft for annotation — *not yet implemented*
**Date:** 2026-08-07
**Target of the documentation:** `https://github.com/cgroza/GraffiTE/tree/v1.1dev` @ `18a76d9`
**Author:** drafted by Claude from a source audit; to be annotated by C. Goubert

> **How to annotate:** sections are numbered so you can reply "§4.2: no, merge these two pages"
> or leave inline `<!-- CG: ... -->` comments. Points that need a decision from you are marked
> **[DECISION]**. Points I have already decided and will just do are marked **[DEFAULT]**.

---

## 1. Why

`README.md` is 935 lines / 84 KB and is the only living documentation. It is thorough, but it is
one flat file that interleaves five different documents:

| Content type | ~lines | % |
|---|---|---|
| Changelog / meta | 170 | 18 % |
| Parameter reference | 215 | 23 % |
| Output & VCF-field reference | 185 | 20 % |
| Biological / method rationale (TSD, HERV-K, L1 inversion, SVA VNTR, `--human`) | 120 | 13 % |
| Tutorial (install, run, update, profiles) | 110 | 12 % |
| Resource tables, FAQ, misc | 135 | 14 % |

Heading levels are inconsistent (`##` jumps to `####` with no `###` under `## Parameters`;
`## TSD module` is a peer of `## Installation`). Two further "docs" exist and are both dead:

- `docs/source/index.rst` — 1066 lines, a Sept-2023 pandoc fork of the README. Three years stale,
  100 % duplicative, wired to no hosting service. Contains `.. raw:: html <details>` blocks that
  survived the conversion. Nothing in it is unique. <!-- Clem: this needs to go... -->
- `README_pipeline_full.md` — a single orphaned Mermaid DAG, unlinked from anything, still showing
  the discontinued `--mammal` branch and missing sniffles2 / PAV / truvari / ULTRA. <!-- Clem: this needs to go... -->

Three goals:

1. A README that sells the tool and routes people, in one screen.
2. A structured, navigable, **verified** reference site. 
3. Output that LLMs (Gemini &c.) can ingest without hallucinating. In practice this means exact
   filenames, explicit defaults with units, source provenance, quarantined design docs, and a
   machine-readable entry point — see §6.

---

## 2. Verification stance

**[DEFAULT]** Nothing is carried over from the current README on trust. Every parameter default,
output filename, INFO tag and process description is re-read from source at `18a76d9`.

An audit pass has already been done (against `f34e921`) and found the README is not merely long but
**wrong in places** — see §9. That is why this is a rewrite against code rather than a
reorganisation of prose.

### 2.1 Blocking prerequisite: sync

The local clone is stale and the working tree is *behind* the documentation target.

```
refs/heads/v1.1dev (actual)   18a76d9   fix(vg_call): ploidy for chrX and chrY
origin/v1.1dev (local ref)    5403488   ← written Aug 7 11:20, predates the PR #92 merge
```

PR **#92** (`v1.1dev-human-filter` → `v1.1dev`, head `f34e921`, +560/−46, 7 files) merged at
18:54 UTC as `dc951bc`. v1.1dev is now three commits beyond that merge:

- `3866f0c` feat: precomputed graph alignments input — adds `params.graph_alignments`
  (CSV `sample,gaf,pack`) and `publishDir "${params.out}/GraffiTE_alignments/"` on `graph_align_reads`
- `32c626b` fix(vg_call): `-a` → `-A` for the giraffe method
- `18a76d9` fix(vg_call): ploidy for chrX and chrY

```bash
git fetch origin
git switch -c docs/restructure origin/v1.1dev
```

The `vg_call`, `graph_align_reads` and precomputed-alignment details in my audit are **known-stale**
and get re-read before anything is written down.

---

## 3. Platform

**[DECIDED — MkDocs Material → GitHub Pages]**

Markdown under `docs/`, built by a GitHub Action, served at `cgroza.github.io/GraffiTE`.

Rationale, and why not the alternatives:

| | MkDocs → Pages | Sphinx → ReadTheDocs | GitHub Wiki |
|---|---|---|---|
| Versioned with the code | ✅ | ✅ | ❌ |
| Docs fix rides the same PR as the code fix | ✅ | ✅ | ❌ |
| Reviewable in PR | ✅ | ✅ | ❌ |
| Source format | Markdown | rST / MyST | Markdown |
| Setup cost | Action + `mkdocs.yml` | RTD account + webhook | none |
| Currently enabled on the repo | n/a | n/a | ❌ (wiki is off) |
| Per-version selector (v1.0 / v1.1) | via `mike` | built-in | ❌ |
| LLM ingestion | plain `.md` in-repo | HTML/rST | HTML only |

The "versioned with the code" row is the decisive one, and §2.1 is a live illustration: the
`--human` filter changed shape entirely in a PR merged today. A wiki would have silently drifted;
in-repo docs fail CI instead.

**[YES]** Do you want the `mike` version selector (so v1.0 docs matching the *Nature Comms*
paper stay live alongside v1.1)? It is ~20 lines of extra config. My inclination: yes, since the
paper is out and people will arrive at v1.0 behaviour. Not in the estimate below if you decline.
<!-- Clem: also preserve the main branch README.md -->
---

## 4. Structure

### 4.1 `README.md` — landing page, target ≤ 150 lines (from 935)

```
Logo + badges  (fix the malformed apptainer badge href, add MIT licence badge)
Paper banner
One-paragraph "what is this"
3 bullets: insertion polymorphisms / VCF annotation / genotyping
Mermaid: three-stage overview
Quickstart: 6 lines, copy-pasteable
Routing table → docs site
Citation · Licence · Issues · Bourque-lab attribution
Star history
```
<!-- Clem: you also need to present the methylation features on the landing page -->

Everything else moves out. The changelog in particular (154 lines, 18 % of the file) goes to
`docs/changelog.md`.

### 4.2 Site tree

```
mkdocs.yml
.github/workflows/docs.yml
docs/
├── index.md                        # what it is, when to use it, the 3-stage model, glossary
├── getting-started/
│   ├── installation.md             # prereqs, --recursive clone, container, profiles
│   ├── quickstart.md               # test dataset, end to end
│   └── choosing-your-inputs.md     # decision tree → which entry flag
├── guides/
│   ├── discovery.md                # Stage A: assemblies / longreads / bams / pav / svs / vcf
│   ├── annotation.md               # Stage B: RepeatMasker + ULTRA + TSD + polyA + filters
│   ├── genotyping.md               # Stage C: pangenie / giraffe / graphaligner / precomputed
│   ├── human-mei.md                # --human pME filter + HERV-K classifier
│   ├── skipping-work.md            # --RM_dir, --graffite_vcf, --graph, --vcfs, --graph_alignments
│   ├── resources.md                # profiles, per-process cpu/mem/time, large genomes
│   └── methylation.md              # --epigenomes / panmethyl
├── reference/
│   ├── parameters.md               # every param: name, default, stage, effect, source line
│   ├── outputs.md                  # full publishDir tree, file by file
│   ├── vcf-fields.md               # every INFO/FORMAT tag + the code that writes it
│   ├── processes.md                # all 21 processes: in / out / tools / resources
│   ├── samplesheets.md             # CSV schema for every --*.csv flag
│   └── container.md                # image contents, pinned vs unpinned
├── background/
│   ├── tsd.md
│   ├── l1-5prime-inversion.md
│   ├── sva-vntr.md
│   └── hervk-hml2.md
├── troubleshooting.md
├── changelog.md
├── design-notes/                   # existing proposals, marked NOT IMPLEMENTED
├── llms.txt
└── assets/
```

**[YES]** Is `background/` the right home for the biology, or would you rather it sat inside
the relevant guide (e.g. TSD rationale inline in `guides/annotation.md`)? Splitting keeps guides
short and is better for retrieval; merging keeps the reader in one place. I lean split, with the
guide linking down.

### 4.3 Files removed or moved

| Action | Path | Why |
|---|---|---|
| **Delete** | `docs/source/`, `docs/Makefile`, `docs/make.bat` | dead Sphinx, nothing unique |
| **Delete** | `README_pipeline_full.md` | superseded by Mermaid diagram §5.3 |
| **Move** | `docs/transduction_module_plan.md`, `..._v2.md`, `docs/cli_ergonomics_and_benchmark_proposal.md`, `docs/HERVK_classification_plan/` | → `docs/design-notes/` |
| **Untouched** | `paper/`, `utils/`, `test/` | content changes only |

**[REPO ONLY]** `docs/design-notes/` — publish them in the site nav (marked "not implemented"), or
keep them in-repo but excluded from the built site? I lean **publish, clearly marked**: they are
good documents and hiding them helps nobody. But see §6 on why the marking matters.

### 4.4 Content sourcing

| Destination | Source of truth |
|---|---|
| `reference/parameters.md` | `nextflow.config` params block + every `params.*` reference in `main.nf` / `module/main.nf` |
| `reference/processes.md`, `outputs.md` | the 21 processes in `module/main.nf`, their `publishDir` and `output:` blocks |
| `reference/vcf-fields.md` | header strings in `bin/repmask_vcf.sh`, `bin/add_polyA.py`, `bin/hervk_classify.py`, `module/main.nf` (`tsd_report`) |
| `reference/container.md` | `GraffiTE.def` |
| `getting-started/choosing-your-inputs.md` | **§1 of `cli_ergonomics_and_benchmark_proposal.md`** — its Stage A/B/C table is a far better mental model than the README's flag list. Reusing the *model*, not the proposals. |
| `background/hervk-hml2.md` | `bin/hervk_classify.py` constants + `HERVK_FILTER_PLAN.md` |
| `reference/vcf-fields.md` (bcftools caveats) | **§1 of `test/human_filter/CLUSTER_TEST.md`** — the `Number=.` element-wise `~` and unreliable-negation discussion is genuine reference knowledge currently trapped in a branch-scoped test doc |
| `background/tsd.md` | `bin/prepTSD.sh`, `bin/TSD_Match_v2.sh`, `bin/exact_match.py` |

<!-- Clem: for HERV-K, trust the code, not the notes that can be outdated -->

---

## 5. Figures and diagrams

**[DECIDED — vendor the PNGs, and add Mermaid]**

All 10 content images are hot-linked from `i.imgur.com` and there is **no local asset directory
anywhere in the repo**. They are unversioned and one imgur policy change from vanishing. They get
copied into `docs/assets/`.

New Mermaid diagrams (text-based; render on GitHub *and* in MkDocs; readable by LLMs as text):

1. **Three-stage overview** — A discovery → B annotation → C genotyping. README + `docs/index.md`.
2. **Entry-point decision tree** — "what data do you have?" → flag → which processes run.
3. **Full process DAG**, all 21 processes. Replaces `README_pipeline_full.md`.
4. **Output tree map** — `out/1_SV_search` … `out/4_Genotyping` + `GraffiTE_graph/` + `GraffiTE_alignments/`.
5. **Filter cascade** — `total_repeat_span` → `pangenome.vcf` → trusted *or* human → HERV-K classes.

Two new hand-drawn SVG cartoons (theme-aware so they work in light and dark):

6. **Reference vs non-reference insertion.** The single biggest source of user confusion: a `DEL`
   record means *TE present in the reference, absent in the sample*. The README already carries a
   prose warning table about TE-presence vs ALT-dosage; a picture fixes it properly.
7. **The `--human` funnel.** `pangenome.vcf` → class gate → subfamily whitelist → SVLEN / ULTRA gate
   → polyA / `n_hits` gate → `pangenome.human.vcf`, with the HERVK+SVA carve-out as a side branch.

**[OK for these two, I will review them]** Any other cartoon you want? Candidates I considered and did not include: the
graph-bubble representation (the existing imgur figure covers it), and a TSD anatomy diagram
(ditto). Say the word and I will add them.

---

## 6. LLM-ingestion measures

Concrete, not vibes:

- **`docs/llms.txt`** — the llms.txt convention: one-line description of the project, then a curated
  link list with a sentence per page.
- **`docs/llms-full.txt`** — every page concatenated into one file, generated by the docs workflow.
  This is the file you paste into Gemini.
- **Provenance on every reference row** — `module/main.nf:247`, `bin/add_polyA.py:MIN_LEN`. A model
  quoting a default can be checked against the line.
- **Explicit defaults with units** — never "default 0.6", always "`0.6` (fraction of variant length)".
- **Version stamp on every page** — *"Applies to GraffiTE v1.1 (`v1.1dev` @ 18a76d9)"*. Prevents a
  model blending v1.0 and v1.1 behaviour, which is a live risk given the paper documents v1.0.
- **Glossary in `docs/index.md`** fixing one term per concept — pME, trusted subset, reference vs
  non-reference insertion, *hit* vs *fragment* — so retrieval does not blend synonyms.
- **Per-page `description` front-matter**; no page depends on a neighbour for its meaning.
- **Design notes quarantined** under `docs/design-notes/` behind an unmissable
  `!!! warning "Not implemented"` admonition. This matters more than it sounds: the transduction
  plans describe INFO fields, new params and a `bin/detect_transductions.py` that **do not exist**.
  Ingested flat, they are exactly what makes a model confidently document a feature you never shipped.

---

## 7. Anti-drift: CI

The README drifted from the code (§9). Docs will too unless something checks. Proposed:

1. **`mkdocs build --strict`** on every PR — fails on any broken internal link or missing nav entry.
2. **`docs/scripts/check_params.py`** — parses the `params { }` block of `nextflow.config` plus every
   `params.X` reference in `main.nf` / `module/main.nf`, diffs against the table in
   `reference/parameters.md`, exits non-zero on mismatch **in either direction**. This is the
   specific defence against findings §9.1 and §9.5 recurring.
3. Same script asserts every `##INFO=<ID=...>` string in `bin/` and `module/main.nf` has a row in
   `reference/vcf-fields.md`.
4. **External-link check** (`lychee`) on a weekly schedule, not per-PR.
   `repeatmasker.org/~cgoubert/GraffiTE_libraries/` is a single-point-of-failure asset worth watching.

**[DECISION]** Is a CI gate welcome, or would you rather these run advisory-only (report, don't
fail)? A hard gate means a param rename cannot merge without a docs update — which is the point,
but it is also friction on a dev branch.

---

## 8. Scope boundary

I will **not** run the pipeline. No claim in the docs will rest on an unverified run. The resource
tables currently in the README (human / *C. sativa* / *Z. mays* CPU-RAM-runtime) will be carried
over **labelled with their provenance and the version they were measured on**, not presented as
current v1.1 figures.
<!-- Clem: yes this table can stay with v1.0 documentation, not in 1.1 -->
---

## 9. Audit findings — for your triage, not silently fixed

Found while verifying. These are code and consistency issues; each needs a call from you. I will
write the docs to describe *what the code does*, and list anything I had to route around.

### 9.1 Documentation is wrong

| # | Finding |
|---|---|
| 1 | **`--pav` is fully wired** (`main.nf:84-88`, `module/main.nf:106`) with **zero** README coverage. Largest single gap. |
| 2 | `mam_filter_1` / `mam_filter_2` documented as live at README:834 and :845, but declared removed at :45-46. Real fields are `L1_5PINV` and a `(VNTR_only)` suffix on `repeat_ids` with the class rewritten to `Simple_repeat`. |
| 3 | README:925 says the TSD module runs only on single-hit SVs; README:708 says otherwise. The code runs it unconditionally. |
| 4 | README documents `--graph_align_theads` (typo; real name `--graph_align_threads`). |
| 5 | Undocumented params: `--tsd_batch_size`, `--tsd_memory`, `--trusted_ignore_filter`, `--pav_memory`, `--pav_time`, `--make_graph_time`, `--vg_call_time`. |

<!-- Clem: OK fix all of that in v1.1 -->

### 9.2 Code / config

| # | Finding |
|---|---|
| 6 | Used in code, never declared in `nextflow.config`: `params.svs`, `params.vcfs`, `params.graph`, `params.bed`, `params.lifted`, `params.graph_alignments`. They work (Groovy null is falsy) but are invisible to `-params-file` and to anyone reading the config. |
| 7 | `main.nf:35` includes `./panmethyl/module/` **unconditionally**, and `panmethyl/` is an uninitialised submodule. **A plain `git clone` cannot run the pipeline.** Docs will say `git clone --recursive`; arguably the include should be conditional. |
| 8 | `GraffiTE.def` does not install `truvari`, `ultra`, or `pyfaidx`, all of which the pipeline calls; `bin/subset_gaf.py` has a `#!/opt/pypy3/bin/pypy3` shebang for a pypy3 that is not installed. The published image and the `.def` have diverged. `reference/container.md` will document the *image*, flagging the `.def` as out of date. |
| 9 | `--vcf` combined with any discovery flag fails: Stage A is skipped when `params.vcf` is set, but `main.nf:109` then reads `sv_variants_ch` → `No such variable`. |
| 10 | `merge_VCFs` uses `publishDir(glob:)`, which is not a publishDir option (it is `pattern:`), and the value `'GraffiTE.merged.genotypes.vcf'` would not match the emitted `.vcf.gz` anyway. |
| 11 | `make_graph`'s graphaligner branch has an unescaped `$PWD` inside a Groovy GString (`export TMPDIR=$PWD`); the giraffe branch correctly uses `\$PWD`. |
| 12 | `graph_method = "precomputed"` is accepted by the branch test in `main.nf` but appears in neither the error message nor the `make_graph` / `graph_align_reads` switches; `nextflow.config`'s comment still lists only three methods. |
| 13 | `bin/TSD_Match_v2.sh` hardcodes the flank width `30` in its scoring awk, so **any `--tsd_win` other than the default silently corrupts the TSD score**. |
| 14 | `pangenie` emits a normalized `${sample_name}.vcf.gz` that its own output glob never captures. |

### 9.3 Repo hygiene

| # | Finding |
|---|---|
| 15 | Dead code in `bin/`: `annotate_vcf_legacy.R`, `one_code_to_find_them_all.pl`, `build_dictionary.pl`, `TSD_Match.sh`, `findTSD.sh`, `findTSD.sh.bak`. All staged into every process by `moduleBinaries`. |
| 16 | Byte-identical duplicates: `HERVK_annot_test.tsv` (`docs/` ↔ `test/`), `repeatmask_gfa.sh` and `join_annotation.R` (`utils/` ↔ `paper/`). |
| 17 | `test/human_test_set.tar.gz` ships a pre-computed `out/` from Sept 2023 that predates ULTRA, `total_repeat_span`, `polyA` and the trusted subsets. **Its "expected output" no longer matches v1.1** — the quickstart cannot honestly tell people to diff against it. |
| 18 | `test/human_filter/run_test.sh` reconstructs the `--human` bcftools expression by `sed`-ing defaults out of `nextflow.config`, deliberately duplicating the Groovy in `concat_repeatmask`. Its own header calls this a sync hazard. |
| 19 | Version is stated three ways: `version.txt` = `1.1.0`, `manifest.version` = `1.1-dev`, `docs/source/conf.py` = `beta`. |
| 20 | No `.github/` at all — no CI, no issue templates, no `CONTRIBUTING.md`. `.DS_Store` files are committed at repo root, `docs/`, `paper/` and three `paper/Applications_Examples/*/`. |

**[DECISION]** For §9.2 and §9.3: do you want me to open GitHub issues, fix them in this same
branch, or leave them entirely and just document around them? My default is **document around
them and hand you the list** — mixing a docs restructure with code fixes makes the PR unreviewable.
Exceptions I would argue for fixing inline because they are one-liners that make the docs honest:
#4 (typo), #6 (declare the missing params), #19 (version consistency).
<!-- Clem: §9.2 no github issue, remove inconsistencies, just write the new doc according to code; v1.0 will inherit whatever is README.md on main branch -->
---

## 10. Proposed order of work

1. Sync to `18a76d9`, branch `docs/restructure` (§2.1). <!-- Clem: no sync to the LATEST commit from the branch v1.1dev whenever you start -->
2. Re-verify the three stale commits; regenerate the process/param/INFO inventories. <!-- Clem: should be fixed -->
3. Scaffold: `mkdocs.yml`, `.github/workflows/docs.yml`, nav, empty pages, `mkdocs build --strict` green.
4. Vendor the imgur assets; write the five Mermaid diagrams and two SVG cartoons.
5. Write `reference/` first — it is the verified backbone everything else links into.
6. Write `guides/`, then `getting-started/`, then `background/`.
7. Rewrite `README.md` last, once the routing targets exist.
8. `docs/llms.txt` + the `llms-full.txt` generator.
9. `check_params.py` + CI wiring.
10. Read the built site end-to-end against `18a76d9`; open the PR. <!-- Clem: against LATEST commit from the branch v1.1dev whenever you start -->

---

## 11. Open questions, collected

- §3 — `mike` version selector for v1.0 vs v1.1 docs? <!-- Clem: ok but what will v1.0 will show? It either need to be built against the main branch, or simply redirecting to the README.md on the current main branch -->
- §4.2 — biology in `background/`, or inline in the guides? <!-- Clem: YES -->
- §4.3 — publish `design-notes/` in the site nav, or repo-only? <!-- Clem: REPO ONLY -->
- §5 — any further cartoons you want? <!-- Clem: YES, see notes -->
- §7 — CI as a hard gate, or advisory? 
- §9 — issues, inline fixes, or hand-off for the audit findings? <!-- Clem: see notes -->
- **Not raised above:** do you want the docs site to carry a "v1.0 (paper) vs v1.1 (current)"
  behaviour-differences page? Readers arriving from the *Nature Comms* paper will hit `--mammal`,
  `mam_filter_*` and OneCode, none of which exist any more. I think this is worth one page. <!-- Clem: YES -->
