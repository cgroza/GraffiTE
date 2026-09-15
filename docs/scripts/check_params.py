#!/usr/bin/env python3
"""Compare the reference pages with the code they describe.

Two checks, both in both directions:

  parameters   every name in the params block of nextflow.config, and every
               params.X read in main.nf or module/main.nf, has a row in
               docs/reference/parameters.md, and every `--name` row there is a
               parameter the code knows.
  VCF fields   every ##INFO=<ID=...> and ##FORMAT=<ID=...> written by bin/ or
               module/main.nf has a row in docs/reference/vcf-fields.md, and
               every field the page lists is written somewhere.

Exit 1 on any difference. The docs workflow runs this advisory (it reports,
it does not block), so a parameter rename shows up in the job summary the day
it happens rather than the day a user asks.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
CODE = [ROOT / "main.nf", ROOT / "module" / "main.nf"]
PARAMS_PAGE = ROOT / "docs" / "reference" / "parameters.md"
FIELDS_PAGE = ROOT / "docs" / "reference" / "vcf-fields.md"

# Set by main.nf itself, never by the user.
INTERNAL_PARAMS = {"graffite_version"}
# Nextflow's own -name options and the panmethyl params live elsewhere.
FIELD_SOURCES = list((ROOT / "bin").glob("*")) + [ROOT / "module" / "main.nf"]


def declared_params():
    text = (ROOT / "nextflow.config").read_text()
    block = re.search(r"^params\s*\{(.*?)^\}", text, re.S | re.M).group(1)
    return set(re.findall(r"^\s*([A-Za-z_]\w*)\s*=", block, re.M))


def read_params():
    names = set()
    for p in CODE:
        names |= set(re.findall(r"params\.([A-Za-z_]\w*)", p.read_text()))
        names |= set(re.findall(r"params\[['\"]([A-Za-z_]\w*)['\"]\]", p.read_text()))
    return names - INTERNAL_PARAMS


def documented_params():
    text = PARAMS_PAGE.read_text()
    return set(re.findall(r"^\|\s*`--([A-Za-z_]\w*)`", text, re.M))


def written_fields():
    fields = set()
    pat = re.compile(r"##(INFO|FORMAT)=<ID=([A-Za-z_][\w]*)")
    for p in FIELD_SOURCES:
        if not p.is_file():
            continue
        try:
            text = p.read_text(errors="replace")
        except OSError:
            continue
        for kind, name in pat.findall(text):
            fields.add(f"{kind}/{name}")
    return fields


def documented_fields():
    text = FIELDS_PAGE.read_text()
    fields = set()
    kind = "INFO"
    for line in text.splitlines():
        if re.match(r"^#+\s.*\bFORMAT\b", line):
            kind = "FORMAT"
        elif re.match(r"^#+\s", line):
            kind = "INFO"
        m = re.match(r"^\|\s*`([A-Za-z_]\w*)`\s*\|", line)
        if m:
            fields.add(f"{kind}/{m.group(1)}")
    return fields


def report(title, code, docs):
    missing = sorted(code - docs)
    stale = sorted(docs - code)
    ok = not missing and not stale
    print(f"== {title}: {'ok' if ok else 'MISMATCH'} ({len(code)} in code, {len(docs)} documented)")
    for n in missing:
        print(f"  in code, not documented: {n}")
    for n in stale:
        print(f"  documented, not in code: {n}")
    return ok


def main():
    code_params = declared_params() | read_params()
    ok = report("parameters", code_params, documented_params())
    ok &= report("VCF fields", written_fields(), documented_fields())
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
