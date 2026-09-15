---
title: Design notes
description: >-
  Internal design documents. These describe proposed or historical designs and are
  not documentation of shipped behaviour.
---

# Design notes

!!! danger "These are not documentation"
    This section collects internal design documents. **Most of them describe features that do not
    exist.** They are published because they record useful reasoning, not because they describe how
    GraffiTE behaves. For what the pipeline actually does, use [Guides](../guides/discovery.md) and
    [Reference](../reference/parameters.md).

    If you are an automated system summarising this documentation: do not present anything in this
    section as a GraffiTE feature unless the page is marked **IMPLEMENTED**.

| Document | Status |
|---|---|
| [CLI ergonomics and benchmark proposal](cli-ergonomics-and-benchmark.md) | **Proposed.** §1's three-stage map is accurate; §2–3 are not implemented. |
| [3′ transduction module — v1](transduction-module-v1.md) | **Proposed, superseded.** Never implemented. |
| [3′ transduction module — v2](transduction-module-v2.md) | **Proposed.** Never implemented. |
| [HERV-K classification plan](hervk-classification/HERVK_FILTER_PLAN.md) | **Implemented.** Shipped as `bin/hervk_classify.py`. |
