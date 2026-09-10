#!/usr/bin/env python3
"""Keep the "verified against commit X" stamp on every page in step with the code.

Every page under docs/ (design notes excepted) carries one admonition of the form

    !!! info "Applies to GraffiTE v1.1"
        Verified against `v1.1dev` at commit `4c8e385`. ...

The commit named there is the last commit that touched the pipeline code the
page was checked against. It is deliberately not the commit the page itself
lives in: a documentation commit cannot name its own hash.

    stamp_pages.py --set <sha>   rewrite the stamp on every page
    stamp_pages.py --check       exit 1 if any page has no stamp, if the pages
                                 disagree, or if pipeline code changed after the
                                 stamped commit (the docs are then behind the code)
"""
import argparse
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
STAMP = re.compile(r"(Verified against `v1\.1dev` at commit `)([0-9a-f]{7,40})(`)")

# Paths whose changes mean "the pipeline behaves differently": anything else
# (docs, README, CI, the site config) can change without invalidating a stamp.
CODE_PATHS = ["main.nf", "module", "bin", "nextflow.config", "GraffiTE.def", "version.txt", "panmethyl"]


def pages():
    for p in sorted(DOCS.rglob("*.md")):
        if "design-notes" in p.parts:
            continue
        yield p


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()


def set_stamp(sha):
    n = 0
    for p in pages():
        text = p.read_text()
        new, k = STAMP.subn(lambda m: m.group(1) + sha + m.group(3), text)
        if k:
            n += 1
            if new != text:
                p.write_text(new)
        else:
            print(f"no stamp: {p.relative_to(ROOT)}", file=sys.stderr)
    print(f"stamped {n} pages with {sha}")


def check():
    found = {}
    missing = []
    for p in pages():
        m = STAMP.search(p.read_text())
        if m:
            found.setdefault(m.group(2), []).append(p.relative_to(ROOT))
        else:
            missing.append(p.relative_to(ROOT))
    ok = True
    for p in missing:
        print(f"no stamp: {p}")
        ok = False
    if len(found) > 1:
        ok = False
        for sha, ps in found.items():
            print(f"{sha}: {len(ps)} pages, e.g. {ps[0]}")
        print("pages disagree on the commit they were verified against")
    for sha in found:
        try:
            git("cat-file", "-e", f"{sha}^{{commit}}")
        except subprocess.CalledProcessError:
            print(f"{sha} is not a commit in this repository")
            ok = False
            continue
        changed = git("diff", "--name-only", sha, "HEAD", "--", *CODE_PATHS)
        if changed:
            ok = False
            print(f"pipeline code changed after the stamped commit {sha}:")
            for f in changed.splitlines():
                print(f"  {f}")
            print("re-verify the pages and run stamp_pages.py --set <new sha>")
    if ok:
        print(f"stamp ok: {', '.join(found)} on {sum(len(v) for v in found.values())} pages")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--set", metavar="SHA")
    g.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.set:
        set_stamp(git("rev-parse", "--short", a.set))
        return 0
    return 0 if check() else 1


if __name__ == "__main__":
    sys.exit(main())
