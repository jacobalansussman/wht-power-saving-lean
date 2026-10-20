#!/usr/bin/env python3
"""d-pair-market: what the U8 word actually does inside ONE Hamming block (default block 0), in the language of the
visit models (pair_toy / m1): which registers sit at which sub-frames of the block span (lines, planes), what they
hold there (restricted to the block's 8 sources), and how targets read block content (pattern, target frame before
the read, whether it is the target's first nonzero frame).
Frames: census climb rule (phases A, centres, precondition, cut, compensation, B, tail, cleanup; reads climb too).
Values: over Q, sources symbolic, helpers/targets start at 0 (dirt is irrelevant for the in-block content), updated
on every add/scale of every phase except centres (cut/compensation only write targets).
Usage: python3 -I u8_block_anatomy.py WORD [block] [out.json]"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
bi = int(sys.argv[2]) if len(sys.argv) > 2 else 0
import os; here = os.path.dirname(os.path.abspath(sys.argv[0]))
blocks = json.load(open(here + '/a1_blocks.json'))['blocks']; S8 = blocks[bi]; S8s = set(S8)
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
span8 = basis([ports[s] for s in S8])
def inblock(f): return all(inspan(x, span8) for x in B[f])
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); T = 0
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
chain = defaultdict(list)
for r in range(n): chain[r].append((-1, dims[cur[r]], inblock(cur[r]) and dims[cur[r]] > 0, None))
def b0(r): return {s: x for s, x in val[r].items() if s in S8s}
def pat(m):
    if not m: return ()
    ss = sorted(m); x0 = m[ss[0]]; return tuple((s, str(m[s]/x0)) for s in ss)
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f; chain[r].append((T, dims[f], inblock(f), pat(b0(r))))
reads = []      # target reads of registers holding nonzero in-block content
firstnz = {}    # target -> time of its first nonzero frame
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        T += 1; continue
    for e in c.get(ph) or []:
        T += 1
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            ta = (dims[cur[a]], inblock(cur[a])); tb = (dims[cur[b]], inblock(cur[b]))
            climb(a, e[-1]); climb(b, e[-1])
            if cls(a) == 'y' and b0(b):
                reads.append(dict(T=T, ph=ph, t=a - v, reg=b, cls=cls(b), pat=pat(b0(b)), other=len(val[b]) - len(b0(b)),
                                  tframe_before=ta, regframe_before=tb, frame_after=dims[cur[a]], first=(ta[0] == 0)))
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den)
            climb(a, e[-1]); val[a] = {key: x*k for key, x in val[a].items()}
for r in range(n): climb(r, c['final'][r])
# registers that sit at an in-block frame of rank 2 or 3
inb = {r: [(d, p) for (t, d, ib, p) in chain[r] if ib and d in (2, 3)] for r in range(n)}
users = {r: x for r, x in inb.items() if x}
print('block', bi, 'sources', S8)
print('registers with in-block rank-2/3 frames, by class and in-block dim path:',
      sorted(Counter((cls(r), tuple(d for d, _ in x)) for r, x in users.items()).items()))
print('in-block content size when arriving at that frame:',
      sorted(Counter((cls(r), d, len(p)) for r, x in users.items() for d, p in x).items()))
print('target reads of in-block content: phase, reg class, |pattern|, other-block terms>0, target frame dim before, first:',
      sorted(Counter((x['ph'], x['cls'], len(x['pat']), x['other'] > 0, x['tframe_before'][0], x['first']) for x in reads).items()))
# clean in-block reads by pattern size, and who reads at frame 0
print('reads at target frame 0 (first reads) with in-block content:', sum(1 for x in reads if x['first']))
# target dim paths
print('target dim paths:', Counter(tuple(d for (_, d, _, _) in chain[t]) for t in range(v, 2*v)).most_common(5))
json.dump(dict(block=bi, users={str(r): x for r, x in users.items()}, reads=reads), open(sys.argv[3] if len(sys.argv) > 3 else here + f'/u8_block{bi}_anatomy.json', 'w'))
