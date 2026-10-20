#!/usr/bin/env python3
"""k8-pairs: per U8 block, what each of its 8 targets' rank-2 stop (the low frame used by LF in pm.py) is spent on:
the partner register's class and in-block content (own-block pair / cross-block pair / other), the pair, and its
sign ratio; plus the target gates y += c y (leaver cells) per block.  Census climb rule, values over Q.
Usage: python3 -I target_stops.py WORD U8_BLOCKS out.json"""
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
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
stop = {}; yy = []; xx = []
def climb(r, f, why):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    if v <= r < 2*v and dims[f] == 2 and r not in stop: stop[r] = why
    cur[r] = f
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f, ('centre',))
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if cls(a) == 'y' and cls(b) == 'y': yy.append((ph, a - v, b - v, str(k), dims[cur[a]], dims[cur[b]]))
            if cls(a) == 'x' and cls(b) == 'x': xx.append((ph, a, b, str(k), dims[cur[a]], dims[cur[b]]))
            cont = {s: x for s, x in val[b].items()}
            ss = sorted(cont)
            why = (ph, cls(b), b, tuple(ss), str(cont[ss[1]]/cont[ss[0]]) if len(ss) == 2 else None, dims[cur[b]])
            climb(a, e[-1], why); climb(b, e[-1], ('partner',))
            if ph not in ('cut', 'compensation'):
                for key, x in list(val[b].items()):
                    y = val[a].get(key, 0) + k*x
                    if y: val[a][key] = y
                    else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); climb(a, e[-1], (ph, 'scale'))
            if ph not in ('cut', 'compensation'): val[a] = {key: x*k for key, x in val[a].items()}
per = defaultdict(Counter); rows = defaultdict(list)
for t in range(v):
    w = stop.get(v + t)
    if w is None: per[blk[t]]['none'] += 1; continue
    ss = w[3]
    kind = 'other'
    if len(ss) == 2:
        kind = 'own-pair' if all(blk[s] == blk[t] for s in ss) else ('cross-pair' if blk[t] in {blk[s] for s in ss} else 'foreign-pair')
    per[blk[t]][(kind, w[1])] += 1
    rows[blk[t]].append(dict(t=t, kind=kind, partner_cls=w[1], partner=w[2], pair=ss, ratio=w[4], phase=w[0]))
print('rank-2 stops of each block\'s 8 targets:')
for b in range(len(blocks)): print('  block', b, dict(per[b]))
yyb = Counter((blk[a], blk[b2]) for ph, a, b2, k, da, db in yy)
print('target gates y += c y:', len(yy), 'by phase', Counter(x[0] for x in yy), 'by (block a, block b):', dict(yyb), 'dims', Counter((x[4], x[5]) for x in yy))
print('source gates x += c x:', len(xx), 'by phase', Counter(x[0] for x in xx), 'blocks', Counter((blk[x[1]], blk[x[2]]) for x in xx))
json.dump(dict(per={b: {str(k): x for k, x in per[b].items()} for b in per}, rows=rows, yy=yy, xx=xx), open(sys.argv[3], 'w'), indent=0)
