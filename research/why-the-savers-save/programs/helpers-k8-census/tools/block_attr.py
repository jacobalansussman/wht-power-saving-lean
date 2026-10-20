#!/usr/bin/env python3
"""k8-pairs: per-block helper attribution of a canonical word on the U8 block partition.
For every helper: the U8 blocks whose sources ever appear in its value (over Q, symbolic sources, all phases except
centres; cut/compensation write targets only), the rank-3 planes it visits (census climb rule) and their block(s),
and its entrance/final dims and width.  Output: per block, helpers that touch only that block (by plane-visiting or
not), and the cross-block helpers.
Usage: python3 -I block_attr.py WORD U8_BLOCKS out.json"""
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
touch = defaultdict(set); planes = defaultdict(list); maxcontent = Counter()
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f
    if dims[f] == 3:
        S = tuple(s for s in range(v) if inspan(ports[s], B[f])); planes[r].append(S)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den); climb(a, e[-1]); climb(b, e[-1])
            if ph in ('cut', 'compensation'): continue
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
            if a >= 2*v:
                touch[a] |= {blk[s] for s in val[a]}; maxcontent[a] = max(maxcontent[a], len(val[a]))
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); climb(a, e[-1])
            if ph not in ('cut', 'compensation'): val[a] = {key: x*k for key, x in val[a].items()}
for r in range(n): climb(r, c['final'][r])
H = range(2*v, n)
per = defaultdict(Counter); cross = Counter(); none = 0
for r in H:
    bs = touch[r] | {blk[s] for S in planes[r] for s in S}
    kind = ('plane' if planes[r] else 'nonplane')
    if len(bs) == 1: per[min(bs)][kind] += 1
    elif not bs: none += 1
    else: cross[(kind, tuple(sorted(bs)), len(bs))] += 1
out = dict(per_block={b: dict(per[b]) for b in range(len(blocks))}, cross={str(k): x for k, x in cross.items()}, no_source=none,
           helpers=R, touch_hist=dict(Counter(len(touch[r]) for r in H)))
print('per block (helpers touching only that block):', {b: dict(per[b]) for b in range(len(blocks))})
print('cross-block helpers:', len(list(cross.elements())), Counter((k[0], k[2]) for k in cross.elements()))
print('no source content:', none, ' touch hist:', out['touch_hist'])
json.dump(out, open(sys.argv[3], 'w'), indent=1)
