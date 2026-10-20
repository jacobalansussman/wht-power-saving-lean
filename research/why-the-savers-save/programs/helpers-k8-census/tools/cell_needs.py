#!/usr/bin/env python3
"""k8-pairs: the pm.py needs (t, block b, |demand_t & b| = 4) that a word serves OUTSIDE pm's in-block row classes:
  cross: some read by t carries exactly 2 sources of b and other-block sources in a row of <= 4 sources (a leaver-cell
         cross row at a cross plane);
  yy:    a FACE need whose b-content reaches t partly through target-target gates (y += y) with 2 of the 4 sources and
         a 2-source pair read (block 0's cell-served faces).
Writes {block: [[t, kind], ...]} for pmk.py's DROP, and per-block counts.
Usage: python3 -I cell_needs.py dlv-style WORD U8_BLOCKS DEMAND out.json"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
d = json.load(open(sys.argv[3])); dem = {int(t): {int(s) for s in dd} for t, dd in d['demand'].items()}
v = c['v']; R = c['R']; n = 2*v + R
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
rd = defaultdict(list)
for ph in ['A', 'precondition', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if v <= a < 2*v and val[b]:
                fp = Counter(blk[s] for s in val[b])
                for bb, x in fp.items(): rd[(a - v, bb)].append((cls(b), x, len(val[b])))
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
out = defaultdict(list); cnt = Counter()
for t in range(v):
    for b in range(len(blocks)):
        S = dem[t] & set(blocks[b])
        if len(S) != 4: continue
        kind = 'own' if blk[t] == b else 'face'
        rs = rd[(t, b)]
        cross = any(x == 2 and 2 < tot <= 4 for cl, x, tot in rs)
        yy = kind == 'face' and any(cl == 'y' and x == 2 for cl, x, tot in rs) and any(cl != 'y' and x == 2 and tot == 2 for cl, x, tot in rs)
        if cross or yy:
            out[b].append([t, kind, 'cross' if cross else 'yy']); cnt[(b, kind, 'cross' if cross else 'yy')] += 1
print('cell-served needs (block, need kind, how):', dict(sorted(cnt.items())))
print('per block:', {b: len(x) for b, x in sorted(out.items())})
json.dump({str(b): x for b, x in out.items()}, open(sys.argv[4], 'w'), indent=0)
