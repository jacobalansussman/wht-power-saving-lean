#!/usr/bin/env python3
"""k8-pairs: the non-plane ('M-copy') helpers: for every helper that visits no rank-3 frame (census climb rule),
its content (sources only, clean replay, cut/compensation skipped) at each B/tail read by a target: block footprint,
and whether the content lies in the reading target's own U8 block.  Summary per block: M-copy helpers whose first
target-read content is inside one block (= that block's own-block second halves).
Usage: python3 -I mcopy.py WORD U8_BLOCKS out.json"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cur = list(c['start']); val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
plane = Counter()
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f
    if dims[f] == 3: plane[r] += 1
reads = defaultdict(list)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            fa = dims[cur[a]]
            climb(a, e[-1]); climb(b, e[-1])
            if ph in ('cut', 'compensation'): continue
            if v <= a < 2*v and b >= 2*v and val[b]:
                reads[b].append((ph, a - v, tuple(sorted(Counter(blk[s] for s in val[b]).items())), fa, len(val[b])))
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); climb(a, e[-1])
            if ph not in ('cut', 'compensation'): val[a] = {key: x*k for key, x in val[a].items()}
for r in range(n): climb(r, c['final'][r])
NP = [r for r in range(2*v, n) if not plane[r]]
per = defaultdict(Counter); shape = Counter(); nread = Counter()
for r in NP:
    rs = reads[r]; nread[len(rs)] += 1
    if not rs: per['unread']['-'] += 1; continue
    ph, t, fp, fa, sz = rs[0]
    own = len(fp) == 1 and fp[0][0] == blk[t]
    shape[(len(fp), sz, own, dims[c['start'][r]])] += 1
    if len(fp) == 1: per[fp[0][0]][('own' if own else 'other-target', sz)] += 1
    else: per['multi'][len(fp)] += 1
print('non-plane helpers:', len(NP), ' reads per helper:', dict(nread))
print('first target-read content (#blocks, #sources, inside reader block, helper entrance dim):', dict(shape))
for b in sorted(per, key=str): print('  block', b, dict(per[b]))
json.dump(dict(per={str(b): {str(k): x for k, x in per[b].items()} for b in per}, shape={str(k): x for k, x in shape.items()},
               reads={r: reads[r] for r in NP}), open(sys.argv[3], 'w'), indent=0)
