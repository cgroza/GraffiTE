#!/usr/bin/env python3
"""Move the file:line citations in the docs from one commit to another.

Every reference row carries a citation such as `module/main.nf:482-483`
or `bin/add_polyA.py:21-23,100-118`, read at the commit the page is stamped
with. When the code moves, the numbers go stale. This script follows each
cited line from the old commit to the new one through the unchanged lines
around it, rewrites the citation when the cited lines themselves are
unchanged, and lists the citations whose lines were edited or deleted,
which is where a human has to look again.

    remap_citations.py <old-sha> <new-sha> [--write]

Without --write it only reports. Pages under design-notes/ are skipped.
"""
import difflib
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
CITE = re.compile(r"`((?:[\w./-]+/)?[\w.-]+\.(?:nf|py|sh|R|config|def|txt|yml|json)):(\d[\d,\s-]*)`")


def git_show(sha, path):
    r = subprocess.run(["git", "show", f"{sha}:{path}"], cwd=ROOT, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else None


def line_map(old, new):
    """old line number (1-based) -> new line number for lines that survived unchanged, else None."""
    m = {}
    sm = difflib.SequenceMatcher(a=old, b=new, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(i2 - i1):
                m[i1 + k + 1] = j1 + k + 1
    return m


def follow(n, lm):
    """Where line n went. An edited line is placed by its nearest unchanged
    neighbours when nothing was inserted or removed between them."""
    if n in lm:
        return lm[n], True
    above = max((k for k in lm if k < n), default=None)
    below = min((k for k in lm if k > n), default=None)
    if above is None or below is None:
        return None, False
    if lm[below] - lm[above] == below - above:
        return n + lm[above] - above, False
    return None, False


def remap_spec(spec, lm):
    """'12-15,20' -> ('14-17,22', ok, edited). ok is False if a cited line was
    removed or sits in a region with insertions; edited is True if a cited
    line survived only by position, with its text changed."""
    parts, ok, edited = [], True, False
    for piece in re.split(r",\s*", spec.strip()):
        if "-" in piece:
            a, b = (int(x) for x in piece.split("-"))
            (na, ea), (nb, eb) = follow(a, lm), follow(b, lm)
            edited |= not (ea and eb)
            if na is None or nb is None or nb - na != b - a:
                ok = False
                parts.append(piece)
            else:
                parts.append(f"{na}-{nb}")
        else:
            n = int(piece)
            nn, e = follow(n, lm)
            edited |= not e
            if nn is None:
                ok = False
                parts.append(piece)
            else:
                parts.append(str(nn))
    return ",".join(parts), ok, edited


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    old_sha, new_sha = sys.argv[1], sys.argv[2]
    write = "--write" in sys.argv
    maps, missing = {}, set()
    moved = changed = same = 0
    stale, review = [], []
    for page in sorted(DOCS.rglob("*.md")):
        if "design-notes" in page.parts:
            continue
        text = page.read_text()

        def sub(m):
            nonlocal moved, changed, same
            path, spec = m.group(1), m.group(2)
            if path not in maps:
                o, n = git_show(old_sha, path), git_show(new_sha, path)
                if o is None or n is None:
                    missing.add(path)
                    maps[path] = None
                else:
                    maps[path] = line_map(o, n) if o != n else "same"
            lm = maps[path]
            if lm is None:
                return m.group(0)
            if lm == "same":
                same += 1
                return m.group(0)
            new_spec, ok, edited = remap_spec(spec, lm)
            if not ok:
                changed += 1
                stale.append(f"{page.relative_to(ROOT)}: {path}:{spec} (lines removed, or inserted around)")
                return m.group(0)
            if edited:
                review.append(f"{page.relative_to(ROOT)}: {path}:{spec} -> {new_spec} (text of a cited line changed)")
            if new_spec != spec:
                moved += 1
            return f"`{path}:{new_spec}`"

        new_text = CITE.sub(sub, text)
        if write and new_text != text:
            page.write_text(new_text)
    print(f"{moved} citations moved, {same} in unchanged files, {len(review)} moved with edited text, {changed} need a human:")
    for s in stale:
        print("  " + s)
    if review:
        print("moved by position; read the new lines once:")
    for s in review:
        print("  " + s)
    for p in sorted(missing):
        print(f"  not found at one of the commits: {p}")
    if not write and moved:
        print("dry run; pass --write to rewrite the moved citations")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
