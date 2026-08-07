import gzip
import sys
import colorsys
import matplotlib.cm as cm

def parse_gfa(path: str):
    """
    Parse a GFA 1 file.

    Returns
    -------
    headers  : list[str]    – raw H lines
    segments : OrderedDict  – id -> {'seq': str, 'tags': list[str]}
    links    : list[dict]   – L records
    paths    : list[dict]   – P records
    other    : list[str]    – every other record type (passed through)
    """
    segments = set()

    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if cols[0] == 'S':
                if len(cols) < 3:
                    sys.exit(f"ERROR line {lineno}: malformed S record: {line!r}")
                segments.add(cols[1])
    return segments


def scalar_to_mono_rgb(
    value: float,
    hue: float = 0.6,
    *,
    as_float: bool = False,
) -> tuple:
    """
    Map a scalar in [0.0, 1.0] to a monochromatic RGB colour.

    Parameters
    ----------
    value    : float  Input scalar, clamped to [0.0, 1.0].
    hue      : float  HSL hue in [0.0, 1.0].
                        0.0 / 1.0 = red   0.33 = green
                        0.6       = blue  0.75 = magenta
    as_float : bool   Return floats in [0.0, 1.0] instead of
                      integers in [0, 255].

    Returns
    -------
    tuple  (r, g, b)
    """
    value = 1.0 - max(0.0, min(1.0, float(value)))

    r, g, b, _alpha = cm.hot(value)

    return f"#{round(r * 255):02X}{round(g * 255):02X}{round(b * 255):02X}"


segments = parse_gfa(sys.argv[1])

cpg_scores = dict()
cpg_n = dict()
cpg_depth = dict()


for path in sys.argv[2:]:
    levels = gzip.open(path, mode = 'rt', encoding='ascii')
    levels.readline()
    for line in levels:
        node, pos, strand, depth, score, graph = line.rstrip().split('\t')
        score = float(score)
        depth = float(depth)
        new_node = node + "_" + pos
        if new_node in segments:
            if new_node in cpg_scores:
                cpg_scores[new_node] = cpg_scores[new_node] + score
                cpg_n[new_node] = cpg_n[new_node] + 1
                cpg_depth[new_node] = cpg_n[new_node] + 1
            else:
                cpg_scores[new_node] = score
                cpg_n[new_node] = 1
                cpg_depth[new_node] = depth

print('Name,Color,PML,PMD')
for node in cpg_scores:
    score = cpg_scores[node]/cpg_n[node]
    depth = cpg_depth[node]
    rgb_score = scalar_to_mono_rgb(score)
    print("{node},{rgb_score},{score},{depth}".format(node=node, pos=pos, rgb_score=rgb_score, score=score, depth=depth))
