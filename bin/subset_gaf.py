#!/opt/pypy3/bin/pypy3
import sys

for line in sys.stdin:
    cs = None
    fields = line.split('\t')

    for f in fields[12:]:
        if f.startswith('cs:Z:') or f.startswith('cg:Z:'):
            cs = f
            break

    # temporary work around (vg convert GAM to GAF sets MAPQ=0)
    basic_fields ='\t'.join(fields[:12])
    if cs is not None:
        print(basic_fields + "\t" + cs)
    else:
        print(basic_fields)
