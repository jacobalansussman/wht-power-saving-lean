#!/usr/bin/env python3
"""v-k8-0: replay the canonical word over Q (phases precondition, A, B, cleanup, as a0_demand.py); for every
target read of a helper that never visits an in-block line/plane (the 112 'M copies'), record which blocks the
read content's sources lie in, relative to the reader target's own block (target t's own block = block holding
label t).  Also: total own-block content delivered to targets by target<->target gates.
Usage: python3 -I mcopy_reads.py WORD u8_blocks.json"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']
blk = {s: i for i, b in enumerate(blocks) for s in b}
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def red(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x
def basis(b):
    bb = []
    for y in b:
        y = red(y, bb)
        if y: bb.append(y); bb.sort(reverse=True)
    return bb
B = [basis(f) for f in F]; dims = [len(b) for b in B]
spans = [basis([ports[s] for s in x]) for x in blocks]
def inblock(f): return dims[f] in (2, 3) and any(all(red(x, sp) == 0 for x in B[f]) for sp in spans)
at = defaultdict(set)
for r in range(n): at[r].add(c['start'][r]); at[r].add(c['final'][r])
for ph in ['A', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        if e[-1] is None: continue
        at[e[1]].add(e[-1])
        if e[0] == 'add': at[e[2]].add(e[-1])
M = {r for r in range(2*v, n) if not any(inblock(f) for f in at[r])}
print('M copies', len(M))
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
stat = Counter(); readers = defaultdict(set); ownM = 0
for ph in ['precondition', 'A', 'B', 'cleanup']:
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if v <= a < 2*v and b in M:
                t = a - v; srcs = [s for s in val[b]]
                bl = Counter(blk[s] for s in srcs)
                stat[(ph, len(srcs), len(bl), blk[t] in bl)] += 1; readers[b].add(t)
                if blk[t] in bl: ownM += 1
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
print('M-copy reads by (phase, #sources, #blocks, reader own block among them):', dict(stat))
print('M copies by number of distinct reader targets:', Counter(len(x) for x in readers.values()))
print('M-copy reads carrying reader own-block sources:', ownM)
