#!/usr/bin/env python3
"""k8-pairs: per U8 block, its 8 targets' first nonzero frame: lifted to rank 1 first (by what partner class and
entrance dim, content size) or straight to rank 2 (first read = the rank-2 stop).  Census climb rule.
Usage: python3 -I first_reads.py WORD U8_BLOCKS"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
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
cur = list(c['start']); first = {}
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres': continue
    for e in c.get(ph) or []:
        regs = [e[1]] + ([e[2]] if e[0] == 'add' else [])
        for r in regs:
            f = e[-1]
            if f is None or cur[r] == f or dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): continue
            if v <= r < 2*v and r not in first and dims[cur[r]] == 0:
                p = e[2] if e[0] == 'add' and r == e[1] else None
                first[r] = (ph, dims[f], cls(p) if p is not None else 'scale/partner', dims[c['start'][p]] if p is not None and p >= 2*v else None)
            cur[r] = f
per = defaultdict(Counter)
for t in range(v): per[blk[t]][first.get(v + t)] += 1
for b in range(len(blocks)): print(' block', b, dict(per[b]))
