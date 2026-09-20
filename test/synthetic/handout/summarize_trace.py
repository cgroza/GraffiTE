#!/usr/bin/env python3
"""Wall time and peak memory per process, from every cell's nextflow_trace.txt.

  summarize_trace.py <WORKDIR>/runs
"""
import glob
import os
import re
import sys


def to_bytes(s):
    s = (s or '').strip()
    m = re.match(r'^([\d.]+)\s*([KMGT]?B?)$', s, re.I)
    if not m:
        return 0
    v = float(m.group(1))
    return int(v * {'': 1, 'B': 1, 'KB': 1 << 10, 'MB': 1 << 20,
                    'GB': 1 << 30, 'TB': 1 << 40}.get(m.group(2).upper(), 1))


def to_ms(s):
    s = (s or '').strip()
    if s in ('-', ''):
        return 0
    total = 0.0
    for v, u in re.findall(r'([\d.]+)\s*(ms|s|m|h|d)', s):
        total += float(v) * {'ms': .001, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}[u]
    return total


def human_t(sec):
    if sec < 60:
        return '%.1fs' % sec
    if sec < 3600:
        return '%dm%02ds' % (sec // 60, sec % 60)
    return '%dh%02dm' % (sec // 3600, (sec % 3600) // 60)


def human_b(b):
    for u in ['B', 'KB', 'MB', 'GB', 'TB']:
        if b < 1024:
            return '%.1f%s' % (b, u)
        b /= 1024.0
    return '%.1fPB' % b


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    rows = []
    for t in sorted(glob.glob(os.path.join(root, '*', 'nextflow_trace.txt'))):
        cell = os.path.basename(os.path.dirname(t))
        with open(t) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                v = line.rstrip('\n').split('\t')
                if len(v) < len(hdr):
                    continue
                d = dict(zip(hdr, v))
                rows.append((cell, d))

    print('%-20s %-26s %-9s %-10s %-10s %-10s %s'
          % ('cell', 'process', 'status', 'realtime', 'duration', 'peak_rss', 'exit'))
    print('-' * 100)
    agg = {}
    for cell, d in rows:
        name = d.get('name', '?')
        proc = name.split(' ')[0].split('(')[0].strip()
        rt = to_ms(d.get('realtime', ''))
        du = to_ms(d.get('duration', ''))
        rss = to_bytes(d.get('peak_rss', ''))
        print('%-20s %-26s %-9s %-10s %-10s %-10s %s'
              % (cell, proc[:26], d.get('status', '?')[:9], human_t(rt), human_t(du),
                 human_b(rss), d.get('exit', '?')))
        k = (cell, proc)
        a = agg.setdefault(k, [0, 0, 0])
        a[0] += 1
        a[1] = max(a[1], rt)
        a[2] = max(a[2], rss)

    print()
    print('== per process, max over tasks ==')
    print('%-20s %-26s %-6s %-10s %s' % ('cell', 'process', 'tasks', 'max_time', 'max_rss'))
    print('-' * 80)
    for (cell, proc), (n, rt, rss) in sorted(agg.items()):
        print('%-20s %-26s %-6d %-10s %s' % (cell, proc[:26], n, human_t(rt), human_b(rss)))


if __name__ == '__main__':
    main()
