#!/usr/bin/env python3
"""g-break: per-plane anatomy of the U8 word (which registers visit each packet plane, who reads them, the
intersection of the readers' u_t^perp).  Usage: python3 -I explore.py WORD"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
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
dot = lambda a, b: bin(a & b).count('1') & 1
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); visits = defaultdict(list); plane_of = {}
step = 0; first_at = {}
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f; visits[r].append((step, f))
reads = defaultdict(list)   # helper -> list of (step, reader, phase)
writes = defaultdict(list)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        step += 1
        climb(e[1], e[-1])
        if e[0] == 'add':
            climb(e[2], e[-1]); reads[e[2]].append((step, e[1], ph)); writes[e[1]].append((step, e[2], ph))
for r in range(n): climb(r, c['final'][r])
P_of = {}
for r in range(n):
    for s, f in visits[r]:
        if dims[f] == 3: P_of[r] = tuple(B[f])
for r in range(n):
    if dims[c['start'][r]] == 3: P_of[r] = tuple(B[c['start'][r]])
planes = sorted(set(P_of.values()))
vis = defaultdict(list)
for r, P in P_of.items(): vis[P].append(r)
# block structure
srcs_in = {P: [s for s in range(v) if inspan(ports[s], list(P))] for P in planes}
blk = defaultdict(list)
for P in planes: blk[tuple(sorted(set(sum((srcs_in[Q] for Q in planes if set(srcs_in[Q]) & set(srcs_in[P])), []))))].append(P)
print('labels of sources of first plane block:')
P0 = planes[0]; S0 = srcs_in[P0]
# find block of P0 (component)
comp = set(S0); ch = True
while ch:
    ch = False
    for P in planes:
        if set(srcs_in[P]) & comp and not set(srcs_in[P]) <= comp: comp |= set(srcs_in[P]); ch = True
comp = sorted(comp); print(comp, [bin(ports[s]) for s in comp], 'rank', len(basis([ports[s] for s in comp])))
bp = [P for P in planes if set(srcs_in[P]) <= set(comp)]
print('planes in block', len(bp))
for P in bp:
    regs = vis[P]
    print('plane', [bin(x) for x in P], 'srcs', srcs_in[P], 'visitors', [(r, cls(r), dims[c['start'][r]], dims[c['final'][r]]) for r in regs])
    for r in regs:
        if cls(r) == 'h':
            rd = Counter((cls(a), ph) for s, a, ph in reads[r]); wr = Counter((cls(b), ph) for s, b, ph in writes[r])
            tg = sorted({a for s, a, ph in reads[r] if cls(a) == 'y'})
            print('   h', r, 'read by', dict(rd), ' writes from', dict(wr), ' target readers', tg)
            # intersection of readers' perps
            if tg:
                us = basis([ports[t - v] for t in tg]); print('      rank of reader labels', len(us), ' -> cap dim', 9 - len(us))
            print('      frame dims path', [dims[f] for s, f in visits[r]])
