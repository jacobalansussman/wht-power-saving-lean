#!/usr/bin/env python3
"""d-pair-market: what each U8 target spends its low frames on: the gate that lifts it 0 -> rank 1 (its first read)
and the gate that lifts it to rank 2 (its rank-2 stop), with the partner register's class, entrance dim, content
(source count, and whether the content lies in the target's own U8 block), and the frame ranks.
Frames: census climb rule; values over Q, symbolic sources, dirt ignored; cut/compensation (dirt/decoder reads at
frame 0) are replayed for frames only.
Usage: python3 -I u8_target_low.py WORD"""
import sys, json, gzip, os
from fractions import Fraction as Fr
from collections import Counter, defaultdict
here = os.path.dirname(os.path.abspath(sys.argv[0]))
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
U8 = json.load(open(here + '/u8_blocks.json'))['blocks']; blk = {s: i for i, x in enumerate(U8) for s in x}
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
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
low = defaultdict(dict)
def climb(r, f, why):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    if v <= r < 2*v and dims[f] <= 3 and dims[f] not in low[r]: low[r][dims[f]] = why
    cur[r] = f
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f, ('centre',))
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            cont = val[b]; t = a - v
            why = (ph, cls(b), dims[c['start'][b]] if cls(b) == 'h' else None, len(cont),
                   (all(blk[s] == blk[t] for s in cont) if cont and cls(a) == 'y' else None), dims[cur[b]])
            climb(a, e[-1], why); climb(b, e[-1], ('partner',))
            if ph not in ('cut', 'compensation'):
                for key, x in list(val[b].items()):
                    y = val[a].get(key, 0) + k*x
                    if y: val[a][key] = y
                    else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); climb(a, e[-1], (ph, 'scale'))
            if ph not in ('cut', 'compensation'): val[a] = {key: x*k for key, x in val[a].items()}
for d in (1, 2, 3):
    print(f'targets lifted to rank {d} by (phase, partner class, partner entrance dim, |content|, own-block content, partner frame dim):')
    for k, x in Counter(low[r].get(d, ('never',)) for r in range(v, 2*v)).most_common(8): print('   ', x, k)
